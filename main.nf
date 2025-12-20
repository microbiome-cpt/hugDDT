#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
    println """
    ==================================================
        Microbiome Metagenomics Analysis Pipeline
       human gut diseases diagnostics tool (hugDDT)
    ==================================================

    Usage:
        nextflow run main.nf [options]

    Options:
        --input_ont      Directory containing ONT FASTQ files
        --input_pod      Directory containing raw signal .pod5 files
        --input_short    Directory containing paired Illumina FASTQ files
        --outdir         Output directory [results/]
        --cpus           Number of threads to use [100]
        --mem         Amount of RAW memory to allocate ['500GB']
        --help           Show this message

    Configurable Options: user can redefine all tool-specific parameters from the pipeline run call
        --chopper_opts   Trimming and filtration options for Chopper ['-q 10 -l 1000']
        --racon_round    Number of Racon polish iterations [2] 

    Profiles:
        --hybrid         Performs hybrid long/short reads correction and assembly 
        --expand         Replaces read-based taxonomy profiling module with more 
                         slow and accurate toolset
        --polybin        Enables parallel binning and bin refinement process via 
                         three independent binners (SemiBin2, COMEBin, Metabinner)  
        --slurm          Runs the pipeline as slurm request
 
    Nextflow parameters:
        -resume          Resume previous interrupted analysis
        -with-conda      Enable implemented conda envs 
        -with-docker     Enable process execution in a Docker containers
        -profile         Choose a configuration profile (slurm, conda, docker) [standard]

    Examples:
        nextflow run main.nf --input_ont reads/ont --input_short reads/illumina --output results --hybrid -profile slurm
        nextflow run main.nf --input_ont reads/ont --output report --expand --chopper_opts '-q 9 -l 500'
    """.stripIndent()
}

include { PORECHOPPER; FASTP; DORADO_CORRECT; FMLRC2} from './modules/preprocess.nf' 
include { KRAKENUNIQ; KAIJU; KRAKENTOOLS; CENTRIFUGER} from './modules/read_based.nf'
include { METAFLYE;  MINIMAP2; HYBRIDSPADES; RACON; MEDAKA; NEXTPOLISH; NEXTPOLISH_LONG; DORADO_POLISH } from './modules/assembly.nf'
include { SEMIBIN2; MULTIBIN } from './modules/binning.nf' 
include { GTDBTK; PROKKA; ABRICATE; REPORTING } from './modules/annotation.nf' 
include { METAQUAST; CHECKM2; COVERM } from './modules/evaluation.nf' 

workflow{
	main:
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Input validation
    if (!params.input_ont && !params.input_pod5) {
        error "Please provide ONT FASTQ files path (--input_ont) or raw signal file (--input_pod5)"
    }
    if (params.hybrid && !params.input_short) {
        error "Please provide Illumina FASTQ files path (--input_short)"
    }
    
    
    //ONT data input. At this point the QC routine should be done 
	ont_reads = params.input_ont ? channel.fromPath("${params.input_ont}/*.{fastq,fq,fastq.gz,fq.gz}") : channel.empty()

    // Paired-end Illumina reads: expect input_short to be a directory with both R1 and R2 files
    illumina_reads = params.input_short ?  channel.fromFilePairs("${params.input_short}/${params.pattern_short}",
            flat: true,
            checkIfExists: true
        ) : channel.empty()


    //MODULE 1. PREPROCESS
    choped = PORECHOPPER(ont_reads)
    //correads1 = DORADO_CORRECT(choped)

    choped = choped
        .map { file ->
            def sample = file.name.replaceFirst(/\.trimmed\.fastq(?:\.gz)?$/,'')
            tuple(sample, file)
        }

    if (params.hybrid) {
        short_trimmed = FASTP(illumina_reads)
        correads2 = FMLRC2(short_trimmed.filtered_pe
        .join(choped, by:0)
        .map { sample, short1, short2, ont -> tuple(sample, short1, short2, ont) }
        )
    } 

    if (params.hybrid) {
        reads_ch = correads2
    } else {reads_ch = choped}
	//reads_ch = params.hybrid ? correads2 : choped
    //    .map { file ->
    //        def sample = file.name.replaceFirst(/\.corrected\.trimmed\.fasta(?:\.gz)?$/,'')
    //        tuple(sample, file)
    //    }

//	reads_ch = choped
//        .map { file ->
//            def sample = file.name.replaceFirst(/\.corrected\.fasta(?:\.gz)?$/,'')
//            tuple(sample, file)
//        }

    // MODULE 2. READ-BASED TAXONOMY IDENTIFICATION 
	ku = KRAKENUNIQ(reads_ch)
    //ku_out = KRAKENTOOLS(ku.classifications.join(ku.report, by:0)
    //    .join(reads_ch, by:0)
    //    .map { sample, kraken, report, reads -> tuple(sample, kraken, report, reads) })
	cfg = CENTRIFUGER(reads_ch)
	kaiju = KAIJU(reads_ch)


    // MODULE 3. ASSEMBLY AND POLISHING
    if (params.hybrid) {
    assembly_ch = HYBRIDSPADES(reads_ch.join(illumina_reads, by:0)
        .map { sample, kreads, ireads1, ireads2 -> tuple(sample, kreads, ireads1, ireads2) })
    } else {
	assembly_ch = METAFLYE(reads_ch)
    }

    racon = RACON(assembly_ch.assembly
        .join(reads_ch, by:0)
        .map { sample, assembly, reads -> tuple(sample, assembly, reads) })
    medaka = DORADO_POLISH(racon.join(reads_ch, by:0)
        .map { sample, racon_assembly, reads -> tuple(sample, racon_assembly, reads) })

    if (params.hybrid) {
        nextpolish_input = medaka.join(reads_ch, by:0)
        .join(illumina_reads, by:0)
        .map { sample, assembly, reads, short_1, short_2 -> tuple(sample, assembly, reads, short_1, short_2) }
        consensus = NEXTPOLISH(nextpolish_input)
    } else {
        consensus = medaka
    }
    
    //MODULE 4. BINNING
    prebin_map_input = consensus
        .join(reads_ch, by:0)
        .map { sample, assembly, reads -> tuple(sample, assembly, reads) } 
    prebin_map = MINIMAP2(prebin_map_input)


    if (params.polybin) {
        multibinning = MULTIBIN(consensus.polished.join(prebin_map.bam, prebin_map.sam, by:0)
        .map { sample, assembly, bam, sam -> tuple(sample, assembly, bam, sam) } )
//        bins = multibinning.bins
//        mags = multibinning.mags
    } else {
        semibin_input = consensus.polished
        .join(prebin_map.bam, by:0)
        .map { sample, assembly, bam -> tuple(sample, assembly, bam) } 
        semibin = SEMIBIN2(semibin_input)
//        bins = semibin.bins
//        mags = semibin.mags
    }

    //MOSULE 5. EVALUATION
    metaquast = METAQUAST(consensus)
    checkm = CHECKM2(semibin.bins)
    coverm = COVERM(semibin.bins
        .join(prebin_map.bam, by:0)
        .map { sample, bins, bam -> tuple(sample, bins, bam) })


    //MODULE 6. ANNOTATION
    gtdb = GTDBTK(semibin.bins)
    prokka = PROKKA(semibin.mags)
    abricate = ABRICATE(semibin.mags
        .join(semibin.bins, by:0)
        .map { sample, mags, bins -> tuple(sample, mags, bins) })

    report = REPORTING(checkm.report
    .join(coverm.genome_coverm, by:0)
    .join(gtdb.gtdbtk_summary, by:0)
    .join(abricate.vf_sum, by:0)
    .join(abricate.hit_list, by:0)
    .join(abricate.blast_report, by:0)
    .map { sample, checkm_rep, coverm_rep, gtdb_class, abr_sum, abr_list, blast_rep -> tuple(sample, checkm_rep, coverm_rep, gtdb_class, abr_sum, abr_list, blast_rep)})



}
