#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//Genome reconstruction quality check

process METAQUAST {
    tag { "${sample_name} - assembly quality assessment" }
    publishDir "${params.outdir}/assembly/${sample_name}/quast", mode: 'copy'

    input:
    tuple val(sample_name), path(assembly)

    output:
//    tuple val(sample_name), path("quast_results/*"), emit: results
    tuple val(sample_name), path("report.html"), emit: quast
 
    script:
    """
    metaquast.py \
        ${assembly} -o ./ \
        --threads ${task.cpus} \
        ${params.quast_opts}
    """
}

process CHECKM2 {
    conda '/home/PAK-CSPMZ/igrabarnik/miniforge3/envs/pipe_help/'
    tag { "${sample_name} - genome completness assessment" }
    publishDir "${params.outdir}/assembly/${sample_name}/checkm", mode: 'copy'

    input:
    tuple val(sample_name), path(bins)

    output:
    tuple val(sample_name), path("checkm2/quality_report.tsv"), emit: report

    script:
    """
    checkm2 predict --input ${bins} \
    -o ./checkm2 \
    -x .fa \
    -t ${task.cpus} \
    --database_path ${params.checkm2_index} \
    ${params.checkm2_opts}
    """
}

process COVERM {
    tag { "${sample_name} - genome coverage assessment" }
    publishDir "${params.outdir}/assembly/${sample_name}/cover", mode: 'copy'

    input:
    tuple val(sample_name), path(bins), path(bam)

    output:
    tuple val(sample_name), path("genome_cover.tsv"), emit: genome_coverm

    script:
    """
    coverm genome \
    -b ${bam} \
    -d ${bins} \
    -x fa \
    -m mean relative_abundance covered_fraction \
    -t ${task.cpus} \
    -o genome_cover.tsv
    """
}


