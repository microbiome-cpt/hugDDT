# Metagenomics Nextflow Pipeline
human gut Disease Diagnostics Tool )))

## Overview
A lineary-modular Nextflow pipelnie for processing and analysis of high-coverage metagenomic sequencing data from NGS and TGS platforms, aimed at Microbiome Profiling and Search for Virulence Factors associated with human medical conditions. A modular Nextflow pipeline for metagenomic analysis of Nanopore (R10) and hybrid Nanopore+Illumina data of human gut microbiome communities.


## Modules
- **Preprocess**: Adapter trimming, quality filtering, correction (ONT/Illumina)
- **Read_Based**: Taxonomic classification (KrakenUniq, Centrifuger, Kaiju)
- **Assembly**: Assembly and polishing (metaFlye/hybridSPAdes, Racon -> Dorado_polish -> NextPolish)
- **Binning**: MAGs binning (SemiBin2)
- **Annotation**: Taxonomic annotation by GTDB-Tk, virulence/plasmid search by Abricate+BLAST
- **Evaluation**: Genome quality and abundance assessment (CheckM2, CoverM, metaQUAST)

![pipeline chart](workflow.png)

## Profiles
> actually working and usable
- `hybrid`: Enables short reads processing and hybrid assembly/correction
- `prokka`: Enables overall functional analysis (slow but highly detailed)

> actually not working, but prepared to be
- `read_expand`: Adds (or replace) alternative classifier (MEGAN7 or MetaMaps)
- `polybin`: Parallel triple binning and refinement

> actually doesn't exist yet, but prepared to be
- `annotate_expand`: Enables more detailed annotation with eggNOG/KEGG/antismash
- `viral_search`: Enables additional viral genomes search and identifictaion

### Essential Run Parameters

> Note that the parameters starting with the "--" refer to the manual params that are set, findable and also changable in the config file; built-in nextflow functions marked with the single "-".

- `--input_ont`: Path to the Nanopore FASTQ files directory
- `--input_short`: Path to the Illumina paired-end FASTQ files directory
- `--outdir`: Path to the output directory where main results would be published
- `-with-conda`: Several processes run under different conda env, it wouldn't work without this nextflow built-in flag
- `-resume`: Most important flag, turn it on to resume the unfinished run or after appending new samples to the analyzed cohort

### Not That Essential Run Parameters
- `--databases`: Path to the directory containing all the databases (`/mnt/raid0/Databases/DDTdb` by default, do not change) 
- `-work-dir`: path to the work directory where all calculations results would be saved; do not change it if possible (default is `/mnt/raid0/Toxins/work_nf`)
- `--cpus`: How many of the CPUs to allocate for the run (250 by default)
- All other default parameters may be set in the `nextflow.config`

## How to Run the Pipeline

First of all, install the developed environments (look for the `ultimative_installer.sh` script) and activate "pipe_develop" one by `conda activate pipe_develop`.
Then you can run the pipeline from any directory by throwing the next commands sequence:
```bash
nextflow run <path to the main.nf script> --input_ont <path to the long reads fastq files> --outdir <path to the outdir> -with-conda --cpus <num of cpus to allocate>
```

#### Examples of the Pipeline Run
###### Simple run from the Pipeline Directory 
```bash
nextflow run main.nf --input_ont path/to_the/data/ONT_fastq --outdir results_20.12 -with-conda
```

###### Simple run with additional parameters
```bash
nextflow run main.nf --input_ont path/to_the/data/ONT_fastq --outdir results_20.12 -with-conda --cpus 128 --chopper_opts '-q 9 -l 500 --maxlength 200000' 
```

###### Hybrid type of run
```bash
nextflow run main.nf --input_ont path/to_the/data/ONT_fastq --input_short /path/to_the/data/illumina --outdir results_22.12 --hybrid -with-conda --cpus 128 --racon_rounds 4
```

###### Complicated type of run (resumed, different workdir, different run directory)
```bash
nextflow run /mnt/raid0/Toxins/hugDDT/main.nf --input_ont ~/path/to_the/data/ONT_fastq --outdir somewhere/results_30.12 -with-conda --prokka --polybin --cpus 258 -work-dir path/to_the/different/place -resume
```

##### QC Shortcut
Run QC before main pipeline if you want:
```bash
nextflow run qc.nf --input_ont path/to_the/data/ONT_fastq --input_short /path/to_the/data/illumina --outdir reports/QC -C qc.config 
```

## Output Directory Structure
- `annotation/{sample_id}`: Abricate and GTDB-Tk results folder
- `assembly/{sample_id}`: All assembly-related files folder, such as: flye contigs, every polisher output assembly, metaQUAST, CoverM and CheckM2 reports
- `binning/{sample_id}`: All MAGs (bins) fasta files (in `output_bins` folder) and binnig info files (contigs sorting, coverage and algorithmic scores) 
- `taxonomy/{sample_id}`: Taxonomical profiling results (from KrakenUniq, Centrifuge and Kaiju separatly); one day there will be some *consensual* tax profile based on them
- `trimmed/{sample_id}`: Empty symbolic links to the trimmed and filtered reads
- `reports/{sample_id}`: Folder with all output informative reports from all processes
- After all successfully completed runs there will be generated three run reports in the output directory root: report_*.html, timeline_*.html and trace_*.txt 

## Final Reports
- `genomes_report.csv`: All found MAGs characteristics: Taxonomic identification, bin's completness, quantitative proportions, genome sizes, virulence factors per genome, etc.
- `virulence_list.csv`: Detailed list of all the virulence factors found in the community, annotated by the bin taxonomy, local BLAST results and measures if identity
- `genome_cover.tsv`: CoverM results with Mean Coverage, Relative Abundance and Covered Fraction measures for all MAGs
- `quality_report.tsv`: Results of CheckM2 analysis with "specific" mode on 
- `report.html`: metaQUAST interactive report 
- Bunch of PNG files with binning and taxonomical analysis results visualisation 

## Special Scripts
final_reports.py, blast.py abd databinning.sh are used in the pipeline to process or generate some data. The user might be interested only in two of them:
- `scripts/reprocess_report.sh`: Used to rerun the Abricate+BLAST process and final reports generation based on published data, with no need to rerun the pipeline, if any changes in the Abricate or BLAST of reports generator scripts were made. Checks for the already reprocessed samples and yet not calculted, but in progress. Needs one positional argument that would locate the outdir:
```bash
bash scripts/reprocess_reports.sh ~/pipe_results/outdir1
``` 
- `ultimative_installer.sh`: simple script to install all the conda environments and third-party tools, making the execution of the pipeline possible from under any user. **MAMBA BEING INSTALLED IS HIGHLY RECOMMENDED!**

## Help
To print help from the command line:
```bash
nextflow run main.nf --help

``` 
> if you want to check wether the pipeline is running or not, use the `top` or `htop` commands in terminal to check the server resources usage
