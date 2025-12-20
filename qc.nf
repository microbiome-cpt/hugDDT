#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage() {
    println """
    ===================================================
        Microbiome Metagenomics Analysis Pipeline
       human gut diseases diagnostics tool (hugDDT)
    ===================================================

    A standalone quality control workflow shortcut for metagenomic reads preanalysis.

    Usage:
        nextflow run qc.nf -c qc.config [options]

    Options:
        --input_ont      Directory containing ONT FASTQ files
        --input_pod      Directory containing raw signal .pod5 files
        --input_short    Directory containing paired Illumina FASTQ files
        --outdir         Output directory [results/]
        --cpus           Number of threads to use [100]
        --mem            Amount of RAW memory to allocate ['500GB']
        --help           Show this message

    Profiles:
        --basecalling    Performs Dorado basecalling in raw signal files 
                         (input_pod must be specified)

    Examples:
        nextflow run qc.nf -c qc.config --input_ont reads/ont --input_short reads/illumina
        nextflow run qc.nf -c qc.config --input_pod reads/pod5 --outdir reads/ont_new --basecalling --cpus 200 --mem '600GB'
    """.stripIndent()
}

process NANOPLOT {
    tag { "ONT read quality metrics" }
    publishDir "${params.outdir}/qc/long/nanoplot", mode: 'copy'
    maxForks 1  // Run one instance that processes all files

    input:
    path reads

    output:
    path "*"

    when:
    params.input_ont

    script:
    """
    NanoPlot \
        --fastq ${reads} \
        --outdir . \
        --threads ${task.cpus} \
        --plots hex dot
    """
}

process FASTQC_LONG {
    tag { "long read quality assessment" }
    publishDir "${params.outdir}/qc/long", mode: 'copy'
    maxForks 1

    input:
    path reads

    output:
    path "*_fastqc.{zip,html}"

    when:
    params.input_ont

    script:
    """
    fastqc \
        -o . \
        -t ${task.cpus} \
        ${reads}
    """
}

process FASTQC_SHORT {
    tag { "short read quality assessment" }
    publishDir "${params.outdir}/qc/short", mode: 'copy'
    maxForks 1

    input:
    path reads

    output:
    path "*_fastqc.{zip,html}"

    when:
    params.hybrid

    script:
    """
    fastqc \
        -o . \
        -t ${task.cpus} \
        ${reads}
    """
}

process MULTIQC_LONG {
    tag { "aggregate long read QC reports" }
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path './long/*'

    output:
    path "multiqc_long.html"
//    path "multiqc_long_data"


    script:
    """
    multiqc . \
        --force \
        --filename multiqc_long \
        --title "Long Reads QC Report"
    """
}

process MULTIQC_SHORT {
    tag { "aggregate short read QC reports" }
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path './short/*'

    output:
    path "multiqc_short.html"
//    path "multiqc_short_data"

    when:
    params.hybrid

    script:
    """
    multiqc . \
        --force \
        --filename multiqc_short \
        --title "Short Reads QC Report"
    """
}

workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Input validation
    if (!params.input_ont && !params.input_short) {
        error "Please provide either ONT (--input_ont) or Illumina (--input_short) reads"
    }

    // Create channels for input reads
    //ONT data input. At this point the QC routine should be done 
	ont_reads = params.input_ont ? channel.fromPath("${params.input_ont}/*.{fastq,fq,fastq.gz,fq.gz}") : channel.empty()

    // Paired-end Illumina reads: expect input_short to be a directory with both R1 and R2 files
    illumina_reads = params.input_short ?  channel.fromFilePairs("${params.input_short}/${params.pattern_short}",
            flat: true,
            checkIfExists: true
        ).map { id, reads -> reads } : channel.empty()


    // Long read QC
    if (params.input_ont) {
        nanoplot_results = NANOPLOT(ont_reads)
        fastqc_long = FASTQC_LONG(ont_reads)

        // Aggregate long read QC results
        MULTIQC_LONG(
            nanoplot_results.mix(
                fastqc_long
            ).collect()
        )
    }

    // Short read QC
    if (params.input_short) {
        fastqc_short = FASTQC_SHORT(illumina_reads)

        // Aggregate short read QC results
        MULTIQC_SHORT(fastqc_short.collect())
    }
}