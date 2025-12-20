resdir=$1
sample_list=$(ls ${resdir}/trimmed | sed -e 's/.trimmed.fastq//g')

export BLASTDB=/mnt/raid0/Databases/DDTdb/blastdb

for i in $sample_list; do
    if [ -d "${resdir}/reports/${i}" ] && [ ! -d "${resdir}/export/${i}" ]; then 
        echo "report for ${i} exists, processing..."
        cd $resdir
        mkdir -p export/${i}

        cd annotation/${i}/abricate
        python ~/pipe_dev/scripts/blastpy.py ${resdir}/binning/${i}/semibin2/output_bins

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
        cp $checkm $coverm $gtdb $abr_sum $abr_list $blast $regions $quast export/${i}

        python ~/pipe_dev/scripts/final_reports.py $checkm $coverm $gtdb $abr_sum $abr_list $blast
        mv virulence_list.tsv genomes_report.csv export/${i}
    else 
        echo "report for ${i} does not exist or already has been reprocessed, eat shit..."
    fi
done