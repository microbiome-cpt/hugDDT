#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Assembly and polishing module (sample-wise output structure)

process METAFLYE {
    tag { "${sample_name} - metaFlye assembly" }
    publishDir "${params.outdir}/assembly/${sample_name}/contigs", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("*.flye.fasta"), emit: assembly
    tuple val(sample_name), path("flye_info/assembly_info.txt"), emit: info

    script:
    """
    flye \
        --nano-hq ${reads} \
        --out-dir flye_info \
        --threads 128 \
        --meta \
        ${params.flye_opts}

    cp flye_info/assembly.fasta ${reads.simpleName}.flye.fasta
    """
}

/*
process RAVEN {
    tag { "${sample_name} - Raven assembly" }
    publishDir "${params.outdir}/assembly/${sample_name}/contigs", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("*.raven.fasta"), emit: assembly

    script:
    """
    raven \
        --threads ${task.cpus} \
        ${params.raven_opts} \
        ${reads} > ${reads.simpleName}.raven.fasta
    """
}
*/

process HYBRIDSPADES {
    tag { "${sample_name} - HybridSPAdes" }
    publishDir "${params.outdir}/assembly/${sample_name}/contigs_hybrid", mode: 'copy'

    input:
    tuple val(sample_name), path(reads), path(short_1), path(short_2)

    output:
    tuple val(sample_name), path("*.contigs.fasta"), emit: assembly
    tuple val(sample_name), path("spades_info/*"), emit: info

    when:
    params.hybrid

    script:
    """
    spades.py \
        --nanopore ${reads} \
        -1 ${short_1} \
        -2 ${short_2} \
        -o spades_info \
        -t ${task.cpus} \
        --meta \
        ${params.spades_opts}

    cp spades_info/contigs.fasta ${sample_name}.contigs.fasta
    """
}

process MINIMAP2 {
    tag { "${sample_name} - mapping reads to assembly" }
    publishDir "${params.outdir}/binning/${sample_name}/mapping"

    input:
    tuple val(sample_name), path(assembly), path(reads)

    output:
    tuple val(sample_name), path("*.bam"), emit: bam
    tuple val(sample_name), path("*.bam.bai"), emit: bai
    tuple val(sample_name), path("*.sam"), emit: sam

    script:
    """
    minimap2 -d ${sample_name}.mmi ${assembly}

    # Map reads
    minimap2 -ax map-ont -t ${task.cpus} ${sample_name}.mmi ${reads} > ${sample_name}.sam 
 
    samtools sort -@ ${task.cpus} -o ${sample_name}.bam ${sample_name}.sam

    # Index BAM
    samtools index -@ ${task.cpus} ${sample_name}.bam

    # Calculate depth
    # samtools depth ${sample_name}.bam > ${sample_name}.depth
    # check assembly coverage
    # bash ../eprst.sh
    """
}

process RACON {
    tag { "${sample_name} - Racon polishing ${params.racon_rounds} rounds" }
    publishDir "${params.outdir}/assembly/${sample_name}/racon"

    input:
    tuple val(sample_name), path(assembly), path(reads)

    output:
    tuple val(sample_name), path("*polished.fasta"), emit: out

    script:
    """
    cp ${assembly} work.fa
    
    for i in \$(seq 1 ${params.racon_rounds}); do
        minimap2 -x map-ont -t ${task.cpus} work.fa ${reads} > map.paf
        racon -t ${task.cpus} ${params.racon_opts} ${reads} map.paf work.fa > polished.tmp
        mv polished.tmp work.fa
    done

    mv work.fa ${sample_name}.polished.fasta
    """   
}

process MEDAKA {
    conda '/home/PAK-CSPMZ/igrabarnik/miniforge3/envs/pipe_help/'
    tag { "${sample_name} - Medaka polishing" }
    publishDir "${params.outdir}/assembly/${sample_name}/medaka", mode: 'copy'

    input:
    tuple val(sample_name), path(assembly), path(reads)

    output:
    tuple val(sample_name), path("consensus.fasta"), emit: polished

    script:
    """
    medaka_consensus \
        -i ${reads} \
        -d ${assembly} \
        -o ./ \
        -t ${task.cpus} ${params.medaka_opts}
    """
//        ${params.medaka_opts}
//    """
}

process DORADO_POLISH {
    tag { "${sample_name} - Dorado polishing" }
    publishDir "${params.outdir}/assembly/${sample_name}/dorado", mode: 'copy'

    input:
    tuple val(sample_name), path(assembly), path(reads)

    output:
    tuple val(sample_name), path("polished_assembly.fasta"), emit: polished

    script:
    """
    # Align reads to a reference using dorado aligner, sort and index
    ${projectDir}/dorado-1.2.0-linux-x64/bin/dorado aligner ${assembly} ${reads} | samtools sort --threads ${task.cpus} > aligned_reads.bam
    samtools index aligned_reads.bam -@ ${task.cpus}

    # Call consensus
    ${projectDir}/dorado-1.2.0-linux-x64/bin/dorado polish aligned_reads.bam ${assembly} --bacteria --ignore-read-groups > polished_assembly.fasta
    """
}

process NEXTPOLISH {
    tag { "${sample_name} - NextPolish" }
    publishDir "${params.outdir}/assembly/${sample_name}/nextpolish", mode: 'copy'

    input:
    tuple val(sample_name), path(assembly), path(reads), path(short_1), path(short_2)

    output:
    tuple val(sample_name), path("*nextpolish.fasta"), emit: polished

    when:
    params.hybrid

    script:
    """
    ls ${reads} > lgs.fofn
    ls ${short_1} ${short_2} > sgs.fofn
    echo '''
[General]
job_type = local
job_prefix = ${sample_name}_NP
task = best
rewrite = yes
rerun = 3
parallel_jobs = 10
multithread_jobs = ${task.cpus}
genome = ${assembly}
genome_size = auto
workdir = ./
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 1000 -bwa

[lgs_option]
lgs_fofn = ./lgs.fofn
lgs_options = -min_read_len 500 -max_depth 1000
lgs_minimap2_options = -x map-ont
''' > run.cfg

nextPolish run.cfg
    """
}


//пока не нужно:
process HOMOPOLISH {
    conda '/home/PAK-CSPMZ/igrabarnik/pipe_dev/homopolish/environment.yml'
    tag { "${sample_name} - Homopolish polishing" }
    publishDir "${params.outdir}/assembly/${sample_name}/homopolish", mode: 'copy'

    input:
    tuple val(sample_name), path(assembly)

    output:
    tuple val(sample_name), path("*_homopolish.fasta"), emit: polished

    script:
    """
    python3 /home/PAK-CSPMZ/igrabarnik/pipe_dev/homopolish/homopolish.py polish \
        -a ${assembly} \
        -s ${params.homopolish_sketch} \
        -o ./ \
        -t ${task.cpus} \
        ${params.homopolish_opts} \

    cp ${sample_name}_homopolish.fasta ../${sample_name}.polished.fasta
    """
}

process PILON {
    tag { "${sample_name} - Pilon polishing" }
    publishDir "${params.outdir}/assembly/${sample_name}/pilon", mode: 'copy'


    input:
    tuple val(sample_name), path(assembly)
    tuple path(short_1), path(short_2)

    output:
    tuple val(sample_name), path("*.pilon.fasta"), emit: polished

    when:
    params.hybrid

    script:
    """
    # Map short reads to assembly
    bwa index ${assembly}
    bwa mem -t ${task.cpus} ${assembly} ${short_1} ${short_2} | \
        samtools sort -@ ${task.cpus} -o illumina_mapped.bam
    samtools index illumina_mapped.bam

    # Run Pilon
    pilon \
        --genome ${assembly} \
        --frags illumina_mapped.bam \
        --output ${assembly.simpleName}.pilon \
        --changes \
        --threads ${task.cpus}

    mv ${assembly.simpleName}.pilon.fasta ${assembly.simpleName}.pilon.fasta
    """
}

process NEXTPOLISH_LONG {
    tag { "${sample_name} - NextPolish" }
    publishDir "${params.outdir}/assembly/${sample_name}/nextpolish", mode: 'copy'

    input:
    tuple val(sample_name), path(assembly), path(reads)

    output:
    tuple val(sample_name), path("*nextpolish.fasta"), emit: polished


    script:
    """
    ls ${reads} > lgs.fofn
    echo '''
[General]
job_type = local
job_prefix = ${sample_name}_NP
task = best
rewrite = yes
rerun = 3
parallel_jobs = 10
multithread_jobs = ${task.cpus}
genome = ${assembly}
genome_size = auto
workdir = ./
polish_options = -p {multithread_jobs}

[lgs_option]
lgs_fofn = ./lgs.fofn
lgs_options = -min_read_len 1000 -max_depth 1000
lgs_minimap2_options = -x map-ont
''' > run.cfg

nextPolish run.cfg
    """
}
