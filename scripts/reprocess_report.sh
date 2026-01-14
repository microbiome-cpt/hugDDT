#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <results_dir> [--noblast]"
    exit 1
fi

resdir="$1"
NOBLAST=false
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ "${2:-}" == "--noblast" ]]; then
    NOBLAST=true
fi

sample_list=$(ls ${resdir}/trimmed | sed -e 's/.trimmed.fastq//g')

export BLASTDB=/mnt/raid0/Databases/DDTdb/blastdb

for i in $sample_list; do
    if [ -d "${resdir}/reports/${i}" ] && [ ! -d "${resdir}/export/${i}" ]; then 
        echo "report for ${i} exists, processing..."
        cd $resdir
        
        abr_list=annotation/${i}/abricate/hit_list.tsv
        abr_sum=annotation/${i}/abricate/vf_summary.tsv 
        blast=annotation/${i}/abricate/blast_report.tsv 
        regions=annotation/${i}/abricate/extracted_regions.fasta
        contigs=annotation/${i}/abricate/extracted_contigs.fasta
        gtdb=annotation/${i}/gtdbtk/gtdbtk.bac120.summary.tsv
        checkm=assembly/${i}/checkm/checkm2/quality_report.tsv
        coverm=assembly/${i}/cover/genome_cover.tsv
        quast=assembly/${i}/quast/report.html
        virlist=reports/${i}/virulence_list.tsv
        genreport=reports/${i}/genomes_report.csv
        binplot=reports/${i}/bins_plot.png
        mkdir -p export/${i}

        if [[ "${NOBLAST}" == false ]]; then
            (
                cd "${resdir}/annotation/${i}/abricate"
                python ${SCRIPT_DIR}/blastpy.py \
                    "${resdir}/binning/${i}/semibin2/output_bins"

                cd $resdir
                python ${SCRIPT_DIR}/final_reports.py $checkm $coverm $gtdb $abr_sum $abr_list $blast
                mv virulence_list.tsv genomes_report.csv bins_plot.png export/${i}
            )
        else
            cp $virlist $genreport $binplot export/${i}
        fi

        cp $checkm $coverm $gtdb $abr_sum $abr_list $blast $regions $quast export/${i}

    else 
        echo "report for ${i} does not exist or already has been reprocessed, eat shit..."
    fi
done
