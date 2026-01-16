#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Genome annotatoins and Functional analysis


process GTDBTK {
    tag { "${sample_name} - taxonomic classification with GTDB-Tk" }
    publishDir "${params.outdir}/annotation/${sample_name}/gtdbtk", mode: 'copy'

    input:
    tuple val(sample_name), path(bins)

    output:
    tuple val(sample_name), path ("gtdbtk.bac120.summary.tsv"), emit: gtdbtk_summary

    script:
    """
    export GTDBTK_DATA_PATH=${params.gtdb_index}
    
    gtdbtk classify_wf \
        --genome_dir ${bins} \
        --out_dir ./ \
        --cpus ${task.cpus} \
        --pplacer_cpus ${params.pplacer_cpus} \
        --mash_db ./mash_db.msh \
        ${params.gtdbtk_opts}
    """
}

process PROKKA {
    tag { "${sample_name} - functional annotation" }
    publishDir "${params.outdir}/annotation/${sample_name}/prokka", mode: 'copy'

    input:
    tuple val(sample_name), path(mags)

    output:
    tuple val(sample_name), path("*.tsv")

    when:
    params.prokka
    
    script:
    """
    for i in ${mags}; do
        prokka \
        --outdir ./ \
        --force \
        --prefix \${i}_prokka \
        --cpus ${task.cpus} \
        --metagenome \
        ${params.prokka_opts} \
                \$i
    done
    """
}

process ABRICATE {
    tag { "${sample_name} - virulence and resistance genes search" }
    publishDir "${params.outdir}/annotation/${sample_name}/abricate", mode: 'copy'

    input:
    tuple val(sample_name), path(mags), path(bins)

    output:
    tuple val(sample_name), path("vf_summary.tsv"), emit: vf_sum
    tuple val(sample_name), path("hit_list.tsv"), emit: hit_list
    tuple val(sample_name), path("blast_report.tsv"), emit: blast_report
    tuple val(sample_name), path("blast_out.tsv"), emit: blast_out

    script:
    """
    export BLASTDB=${params.databases}/blastdb
    
    # Search against multiple databases
    for i in ${mags}; do
        abricate --db csp_vf --threads ${task.cpus} \$i > \${i}_abr_hits.tsv
    done
    # Summarize
    abricate --summary *.tsv > vf_summary.tsv

    head -n 1 \$(ls *_abr_hits.tsv | head -n1) > hit_list.tsv
    # Добавить все ряды без шапки
    for f in *_abr_hits.tsv; do
        tail -n +2 "\$f" >> hit_list.tsv
    done

    python ${projectDir}/scripts/blastpy.py ${bins}
    """
}

process REPORTING {
    publishDir "${params.outdir}/reports/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(checkm), path(coverm), path(gtdb), path(abr_sum), path(abr_list), path(blast)

    output:
    tuple val(sample_name), path('virulence_list.tsv'), path('genomes_report.csv'), path('bins_plot.png')

    script:
    """
    python ${projectDir}/scripts/final_reports.py ${checkm} ${coverm} ${gtdb} ${abr_sum} ${abr_list} ${blast}
    bash ${projectDir}/scripts/reprocess_reports.sh ${params.outdir} --noblast
    """
}


