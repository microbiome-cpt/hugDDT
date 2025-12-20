#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// Read-based taxonomic classification module (sample-wise output structure)

process KRAKENUNIQ {
    tag { "${sample_name} - KrakenUniq" }
    publishDir "${params.outdir}/taxonomy/${sample_name}/krakenuniq", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("*.kraken_report.tsv"), emit: report
//    tuple val(sample_name), path("*.kraken"), emit: classifications

    when:
    !params.read_expand

    script:
    """
    krakenuniq \
        --db ${params.ku_index} \
        --threads ${task.cpus} --preload \
        --report-file ${reads.simpleName}.kraken_report.tsv \
        --output ${reads.simpleName}.kraken \
        ${reads}
    """
}

process KRAKENTOOLS {
    tag { "${sample_name} - kraken throughput reads extraction" }
    publishDir "${params.outdir}/taxonomy/${sample_name}/rare_reads", mode: 'copy'

    input:
    tuple val(sample_name), path(kraken), path(report), path(reads)

    output:
    tuple val(sample_name), path("*.rare_reads.fasta.gz"), emit: krakreads

    when:
    !params.read_expand

    script:
    """
    ~/pipe_dev/KrakenTools/extract_kraken_reads.py \
    -k ${kraken} \
    -s ${reads} \
    -o ${sample_name}.rare_reads.fasta \
    -t 0 1 2157 \
    -r ${report} \
    --include-children

    gzip ${sample_name}.rare_reads.fasta
    """
}


process CENTRIFUGER {
    tag { "${sample_name} - Centrifuger" }
    publishDir "${params.outdir}/taxonomy/${sample_name}/centrifuger", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("*.cfg_output.tsv"), emit: report
    tuple val(sample_name), path("*.cfg_kreport"), emit: kreport

    when:
    !params.read_expand

    script:
    """
    centrifuger \
        -x ${params.cfg_index} \
        -u ${reads} \
        -t ${task.cpus} ${params.centrifuge_opts} > ${reads.simpleName}.cfg_output.tsv
    centrifuger-kreport \
        -x ${params.cfg_index} \
        --no-lca ${reads.simpleName}.cfg_output.tsv > ${reads.simpleName}.cfg_kreport
    """
}

process KAIJU {
    tag { "${sample_name} - Kaiju" }
    publishDir "${params.outdir}/taxonomy/${sample_name}/kaiju", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)

    output:
    //tuple val(sample_name), path("*.kaiju.tsv"), emit: report
    //tuple val(sample_name), path("*.kaiju.names"), emit: names
    tuple val(sample_name), path("*_kaiju.summary*.tsv"), emit: summary

    when:
    !params.read_expand

    script:
    """
    kaiju \
        -t ${params.kaiju_index}/nodes.dmp \
        -f ${params.kaiju_index}/kaiju_db_refseq.fmi \
        -i ${reads} \
        -z ${task.cpus} \
        -o ${reads.simpleName}.kaiju.tsv

    kaiju-addTaxonNames \
        -t ${params.kaiju_index}/nodes.dmp \
        -n ${params.kaiju_index}/names.dmp \
        -i ${reads.simpleName}.kaiju.tsv \
        -o ${reads.simpleName}.kaiju.names

    kaiju2table \
        -t ${params.kaiju_index}/nodes.dmp \
        -n ${params.kaiju_index}/names.dmp \
        -r species -c 100 \
        -o ${sample_name}_kaiju.summary_S.tsv \
         ${reads.simpleName}.kaiju.names

    kaiju2table \
        -t ${params.kaiju_index}/nodes.dmp \
        -n ${params.kaiju_index}/names.dmp \
        -r genus -c 100 \
        -o ${sample_name}_kaiju.summary_G.tsv \
         ${reads.simpleName}.kaiju.names
    """
}

process METAMAPS {
    tag { "${sample_name} - MetaMaps" }
    publishDir "${params.outdir}/taxonomy/${sample_name}/metamaps", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("*.EM.WIMP"), emit: report

    when:
    params.read_expand

    script:
    """
    metamaps mapDirectly \
        -r ${params.metamaps_index}/DB.fa \
        -q ${reads} \
        -o ${reads.simpleName}.map \
        --threads ${task.cpus} ${params.metamaps_opts}
    metamaps classify \
        --mappings ${reads.simpleName}.map \
        --DB ${params.metamaps_index} \
        -t ${task.cpus}
    """
}

process DIAMOND {
    tag { "${sample_name} - DIAMOND" }
    publishDir "${params.outdir}/taxonomy/${sample_name}/diamond", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("*.diamond.daa"), emit: daa

    when:
    params.read_expand

    script:
    """
    diamond blastx \
        --db ${params.diamond_index} \
        --query ${reads} \
        --out ${reads.simpleName}.diamond.daa \
        --threads ${task.cpus} \
        --outfmt 100 \
        --sensitive
    """
}

process MEGAN_LR {
    tag { "${sample_name} - MEGAN-LR" }
    publishDir "${params.outdir}/taxonomy/${sample_name}/meganlr", mode: 'copy'

    input:
    tuple val(sample_name), path(daa)

    output:
    tuple val(sample_name), path("*.megan.txt"), emit: classifications
    tuple val(sample_name), path("*.rma6"), emit: rma

    when:
    params.read_expand

    script:
    """
    daa2rma \
        --in ${daa} \
        --out ${daa.simpleName}.rma6 \
        --reads ${daa.simpleName}.megan.txt \
        --mapDB ${params.megan_index}
    """
}

process MERGE_TAXONOMY {
    tag { "${sample_name} - merge taxonomy" }
    publishDir "${params.outdir}/taxonomy/${sample_name}/merged", mode: 'copy'

    input:
    tuple val(sample_name), path(classification_files)

    output:
    tuple val(sample_name), path("merged_taxonomy.tsv"), emit: merged

    script:
    """
    echo "Merged taxonomy" #> merged_taxonomy.tsv
    #cat ${classification_files} >> merged_taxonomy.tsv
    """
}

workflow {
    //if (!params.databases) {
    //    error "Please provide the path to taxonomy databases (--databases)"
    //}

    reads_ch = Channel.fromPath(params.hybrid ? "${params.raw_data}/hybrid/*.fastq.gz" : "${params.raw_data}/trimmed/*.fastq.gz")
        .map { file -> tuple(file.simpleName, file) }

    if (params.read_expand) {
        centrifuger_results = CENTRIFUGER(reads_ch)
        diamond_results = DIAMOND(reads_ch)
        megan_results = MEGAN_LR(diamond_results.daa)
        metamaps_results = METAMAPS(reads_ch)

    } else {
        krakenuniq_results = KRAKENUNIQ(reads_ch)
        kaiju_results = KAIJU(reads_ch)
    }
}
