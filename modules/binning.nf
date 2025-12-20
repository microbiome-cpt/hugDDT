#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Binning module for metagenomics pipeline

process SEMIBIN2 {
    tag { "${sample_name} - binning with SemiBin2" }
    publishDir "${params.outdir}/binning/${sample_name}/semibin2", mode: 'copy'

    input:
    tuple val(sample_name), path(assembly), path(bam)

    output:
    tuple val(sample_name), path("output_bins/*.fa"), emit: mags
    tuple val(sample_name), path("output_bins/"), emit: bins
    tuple val(sample_name), path("output_bins/unbinned.fa"), emit: unbinned
    tuple val(sample_name), path("*.tsv")
    tuple val(sample_name), path("*.csv")


    script:
    """
    export LD_LIBRARY_PATH=''

    SemiBin2 single_easy_bin \
        -i ${assembly} \
        -b ${bam} \
        -o ./ \
        -t ${task.cpus} \
        ${params.semibin_opts}

    awk '\$2 == -1 {print \$1}' contig_bins.tsv > unbinned_contigs.txt
    seqkit grep -f unbinned_contigs.txt ${assembly} > output_bins/unbinned.fa
    """
}

process MULTIBIN {
    tag { "${sample_name} - Parallel expanded binning" }
    publishDir "${params.outdir}/binning/${sample_name}/expand_bins", mode: 'copy'

    input:
    tuple val(sample_name), path(assembly), path(bam), path(sam)

    output:
    tuple val(sample_name), path("resulted_bins/*.fa"), emit: mags
    tuple val(sample_name), path("resulted_bins/"), emit: bins


    when:
    params.polybin

    script:
    """
    bash databinning.sh  \
        -m auto \
        -a ${assembly} \
        -b ${bam} \
        -s ${sam} \
        -o ./ \
        -t ${task.cpus} \
     """
}

process METABINNER {
    tag { "binning with MetaBinner" }
    publishDir "${params.outdir}/binning/metabinner", mode: 'copy'

    input:
    path assembly
    path bam

    output:
    path "metabinner_bins/*.fa", emit: bins

    when:
    params.polybin

    script:
    """
    metabinner.sh \
        -a ${assembly} \
        -o metabinner_bins \
        -t ${task.cpus} \
        ${params.metabinner_opts}
    """
}

process COMEBIN {
    tag { "binning with ComeBin" }
    publishDir "${params.outdir}/binning/comebin", mode: 'copy'

    input:
    path assembly
    path bam

    output:
    path "comebin_bins/*.fa", emit: bins

    when:
    params.polybin

    script:
    """
    comebin.py \
        --contigs ${assembly} \
        --bam ${bam} \
        --output comebin_bins \
        --threads ${task.cpus} \
        ${params.comebin_opts}
    """
}

process METAWRAP_REFINE {
    tag { "refining bins with metaWRAP" }
    publishDir "${params.outdir}/binning/metawrap", mode: 'copy'

    input:
    path '*'

    output:
    path "refined_bins/*.fa", emit: bins
    path "stats/*", emit: stats

    when:
    params.polybin

    script:
    """
    metawrap bin_refinement \
        -o . \
        -A metabinner_bins \
        -B comebin_bins \
        -c 50 \
        -x 10 \
        -t ${task.cpus}
    """
}

