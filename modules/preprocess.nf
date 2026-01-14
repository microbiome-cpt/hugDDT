#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Preprocessing module for metagenomics pipeline

process DORADO {
    tag { "basecalling_and_processing" }
    publishDir "${params.outdir}/basecalled", mode: 'copy'
    
    input:
    path fast5_dir

    output:
    path "*.fastq.gz", emit: fastq
    path "sequencing_summary.txt", optional: true, emit: summary

    when:
    params.input_fast5

    script:
    def demux_cmd = params.demultiplex ? "--kit ${params.barcode_kit} --trim_barcodes" : ""
    def mod_cmd = params.detect_modifications ? "--modified-bases ${params.modification_model}" : ""
    def duplex_cmd = params.duplex ? "--duplex" : ""
    """
    # Basecalling with optional demultiplexing and modification detection
    dorado basecaller ${params.dorado_model} \
        ${demux_cmd} \
        ${mod_cmd} \
        ${duplex_cmd} \
        ${params.dorado_opts} \
        ${fast5_dir} > basecalled.fastq

    # If demultiplexing was performed, process each barcode separately
    if [ "${params.demultiplex}" = "true" ]; then
        mkdir -p barcodes
        cat basecalled.fastq | awk -v FS=" " '{
            if (substr(\$0,1,1)==">") {
                split(\$2,a,"=")
                file="barcodes/barcode_"a[2]".fastq"
            }
            print \$0 > file
        }'
        for f in barcodes/*.fastq; do
            gzip \$f
        done
    else
        gzip basecalled.fastq
    fi
    """
}

//process NANOFILT {
//    tag { "quality filtering" }
//    publishDir "${params.outdir}/filtered", mode: 'copy'
//
//    input:
//    path reads
//
//    output:
//    path "*.filtered.fastq.gz", emit: filtered
//
//    script:
//    """
//    zcat ${reads} | NanoFilt ${params.nanofilt_opts} | gzip > ${reads.simpleName}.filtered.fastq.gz
//    """
//}

process PORECHOPPER {
    tag { "adapter trimming and filtering" }
    publishDir "${params.outdir}/trimmed"

    input:
    path reads

    output:
    path "*.trimmed.fastq", emit: trimmed

    script:
    """
    seqkit split2 -s 1000000 -j ${task.cpus} -O tmp_chuncks ${reads}

    ls tmp_chuncks/*.fastq* \
    | parallel -j 32 \
    'porechop -i {} ${params.porechop_opts} -t 16 | \
    chopper -t ${task.cpus} ${params.chopper_opts} > {.}_tr.fastq'

    cat tmp_chuncks/*_tr.fastq > ${reads.simpleName}.trimmed.fastq
    rm -r tmp_chuncks 
    """
//    gzip > ${reads.simpleName}.trimmed.fastq.gz
//    """
}

process FASTP {
    tag { "illumina preprocessing" }
    publishDir "${params.outdir}/trimmed/illumina"

    input:
    tuple val(sample_name), path(short_1), path(short_2)

    output:
    tuple val(sample_name), path("*R1_filtered.fastq.gz"), path("*R2_filtered.fastq.gz"), emit: filtered_pe
    path "fastp.json", emit: json
    path "fastp.html", emit: html

    when:
    params.hybrid && params.input_short

    script:
    """
    fastp -i ${short_1} -I ${short_2} \
        -o ${sample_name}.R1_filtered.fastq.gz \
        -O ${sample_name}.R2_filtered.fastq.gz \
        -w ${task.cpus} \
        ${params.fastp_opts} \
    """
}

process DORADO_CORRECT {
    tag { "ONT self-correction" }
    publishDir "${params.outdir}/trimmed/corr/ont"

    input:
    path ont_reads

    output:
    path "*.corrected.fasta.gz", emit: corrected
    
    when:
    !params.hybrid 

    script:
    """
    export LD_LIBRARY_PATH=~/pipe_dev/dorado-1.3.0-linux-x64/lib:\$LD_LIBRARY_PATH

    ${projectDir}/dorado-1.2.0-linux-x64/bin/dorado correct \
    -t ${task.cpus} ${ont_reads} | pigz -p 100 > ${ont_reads.simpleName}.corrected.fasta.gz
    """
}

process FMLRC2 {
    tag { "hybrid correction" }
    publishDir "${params.outdir}/trimmed/corr/hybrid"

    input:
    tuple val(sample_name), path(short_1), path(short_2), path(ont_reads)

    output:
    tuple val(sample_name), path("*.hybrid.corrected.fastq.gz"), emit: corrected

    when:
    params.hybrid 

    script:
    """
    gunzip -c ${short_1} ${short_2} | \
    awk 'NR % 4 == 2' | \
    tr NT TN | \
    ropebwt2 -LR | \
    tr NT TN | \
    fmlrc2-convert comp_msbwt.npy


    fmlrc2 -t ${task.cpus} comp_msbwt.npy \
    ${ont_reads} ${sample_name}.hybrid.corrected.fasta
    
    pigz -p ${task.cpus} ${sample_name}.hybrid.corrected.fasta
    """
}

workflow {
    //main:
    // Input validation
    if (!params.input_ont && !params.input_fast5) {
        error "Please provide either ONT FASTQ (--input_ont) or FAST5 (--input_fast5) files"
    }

    // Paired-end Illumina reads: expect input_short to be a directory with both R1 and R2 files
    if (params.input_short) {
        illumina_reads = Channel.fromFilePairs(
            "${params.input_short}/${params.pattern_short}",
            flat: true,
            checkIfExists: true
        )
    } else {
        illumina_reads = Channel.empty()
    }

    // Main workflow
    if (params.input_fast5) {
        ont_reads = DORADO(params.input_fast5)
    } else {
        ont_reads = params.input_ont ? Channel.fromPath("${params.input_ont}/*.{fastq,fq,fastq.gz,fq.gz}") : Channel.empty()
    }

    filtered_reads = PORECHOPPER(ont_reads)
    //filtered_reads = NANOFILT(trimmed_reads)

    if (params.hybrid && params.input_short) {
        illumina_processed = FASTP(illumina_reads)
        corrected_reads = FMLRC2(filtered_reads, illumina_processed.filtered_pe)
    } else {
        corrected_reads = filtered_reads
    }

    //emit:
    //filtered_fastq = filtered_reads
    //corrected_fastq = corrected_reads
    //short_reads = params.hybrid ? illumina_processed.filtered_pe : null
} 
