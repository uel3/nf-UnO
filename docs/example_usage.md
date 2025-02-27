## How to run nf-UnO with test data presented as results in example MultiQC output

nf-UnO was run on an Salmonella outbreak set of data which can be accessioned from the Sequence Read Archive using the following IDs:

SRR5058919
SRR5058920
SRR5058921
SRR5058922
SRR5058923
SRR5058924
SRR5058925
SRR5058926
SRR5058927
SRR5058928

## Preparing a samplesheet for running nf-UnO

The salmonella outbreak set of data consists of 10 sets of paired end (PE) reads. These reads should be passed to the nf-UnO using a samplesheet.csv with 4 requuired columns, the 5th column is not required for this version of the nf-UnO pipeline. Future versions will include added support for long reads. 

The 4 required are: sample,group,short_reads_1,short_reads_2,long_reads
sample should correspond to the PE read in some way the user can idenitfy.
group is an arbitrary label used by nf-UnO to denote related samples for co-assembly and read mapping. Related samples should have the same group for proper analysis. 
short_reads_1 is the path to read one of the PE read set.
short_reads_2 is the path to read two of the PE read set.

The samplesheet.csv should be formatted like so:

sample,group,short_reads_1,short_reads_2,long_reads
SRR5058919,UNO,/path/to/reads/SRR5058919_1.fastq.gz,/path/to/reads/SRR5058919_2.fastq.gz,
SRR5058920,UNO,/path/to/reads/SRR5058920_1.fastq.gz,/path/to/reads/SRR5058920_2.fastq.gz,
SRR5058921,UNO,/path/to/reads/SRR5058921_1.fastq.gz,/path/to/reads/SRR5058921_2.fastq.gz,
SRR5058922,UNO,/path/to/reads/SRR5058922_1.fastq.gz,/path/to/reads/SRR5058922_2.fastq.gz,
SRR5058923,UNO,/path/to/reads/SRR5058923_1.fastq.gz,/path/to/reads/SRR5058923_2.fastq.gz,
SRR5058924,UNO,/path/to/reads/SRR5058924_1.fastq.gz,/path/to/reads/SRR5058924_2.fastq.gz,
SRR5058925,UNO,/path/to/reads/SRR5058925_1.fastq.gz,/path/to/reads/SRR5058925_2.fastq.gz,
SRR5058926,UNO,/path/to/reads/SRR5058926_1.fastq.gz,/path/to/reads/SRR5058926_2.fastq.gz,
SRR5058927,UNO,/path/to/reads/SRR5058927_1.fastq.gz,/path/to/reads/SRR5058927_2.fastq.gz,
SRR5058928,UNO,/path/to/reads/SRR5058928_1.fastq.gz,/path/to/reads/SRR5058928_2.fastq.gz,

## Running nf-UnO
This outbreak set of data was run using the following command 

nextflow run nf-uno -profile conda -c cdc-dev.config --midas2_uhgg_db /path/to/local/midasdb/my_midasdb_uhgg --input /path/to/samplesheet/coalSRR_samplesheet.csv --outdir results_coalSRR

-profile conda is required for running MIDAS2 workflow. All processes have support for conda, singularity, and docker.
-c is where institutional configuration files are passed. This configuration file contains information for nf-uno about system resources. Information in the config used to run thise outbreak set include:
--midas2_uhgg_db is for passing a local copy of the midas2 uhgg database, this is not requried but useful in the case where a local copy exists and new download is not necessary.
--input is the path to the propeorly formatted samplesheet.csv (required)
--outdir is the user named output directory for results (required)

### profile configs
Nextflow provides common Nextflow pipeline configurations for comput environments [https://nf-co.re/configs/]
Some important aspects for configuration files include, resoruce allocaiton for your computing environment including but not limited to hpc configuration, container settings (conda, singlarity, and docker), and memory allocation.
Specifically for container configurations, I load miniconda for the conda package manager. Configuration settings for running -profile conda include:
profiles {
  conda {
    conda.enabled = true
    conda.useMamba = false
    conda.cacheDir =  "${HOME}/my_conda_envs/nextflow"
    conda.createTimeout    = "120 min"
    docker.enabled         = false
    singularity.enabled    = false
    podman.enabled         = false
    shifter.enabled        = false
    charliecloud.enabled   = false}
    }
## MultiQC output 
In the docs folder [https://github.com/uel3/nf-UnO/tree/master/docs] there is a MultiQC Report for nf-uno analysis of the salmonella outbreak data set. 
Below is a breakdown of each section of the MultiQC Report:

### General Stats

### MIDAS2 Species Abundance

### FastQC: Raw Reads

### Bowtie2: Host Removal 

### FastQc: After Preprocessing

### Binning: Sample Coverage of DASTool Bins

### CheckM: Bin Quality Statistics

### nf-uno Workflow Summary

### nf-uno Software Versions

### Software Versions

### nf-uno Methods Description 

