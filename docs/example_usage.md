## How to run nf-UnO with test data presented as results in example MultiQC output

nf-UnO was run on an Salmonella outbreak set of data which can be accessioned from the Sequence Read Archive using the following IDs:

```console
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
```

## Preparing a samplesheet for running nf-UnO

The salmonella outbreak set of data consists of 5 sets of paired end (PE) reads. These reads should be passed to the nf-UnO using a samplesheet.csv with 4 required columns, the 5th column is not required for this version of the nf-UnO pipeline. Future versions will include added support for long reads. 

The 4 required are: 

```console
sample,group,short_reads_1,short_reads_2,long_reads
```
`sample` should correspond to the PE read in some way the user can idenitfy.

`group` is an arbitrary label used by nf-UnO to denote related samples for co-assembly and read mapping. Related samples should have the same group for proper analysis. 

``short_reads_1`` is the path to read one of the PE read set.

``short_reads_2`` is the path to read two of the PE read set.

The samplesheet.csv should be formatted like so:

```console
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
```

## Running nf-UnO
This outbreak set of data was run using the following command: 

```console
nextflow run nf-uno -profile conda -c cdc-dev.config --midas2_uhgg_db /path/to/local/midasdb/my_midasdb_uhgg --input /path/to/samplesheet/colSRR_samplesheet.csv --outdir results_colSRR
```

``-profile`` conda is required for running MIDAS2 workflow. All processes have support for conda, singularity, and docker.

``-c`` is where institutional configuration files are passed. This configuration file contains information for nf-uno about system resources. Information in the config used to run thise outbreak set include:

``--midas2_uhgg_db`` is for passing a local copy of the midas2 uhgg database, this is not required but useful in the case where a local copy exists and new download is not necessary. If using the midas2db downloaded after an initial run of nf-uno, the name is set to my_midasdb_uhgg by the midasdb process. Existing midas2 database directories may have custom names and provided accordingly (optional)

``--input`` is the path to the propeorly formatted samplesheet.csv (required)

``--outdir`` is the user defined output directory for results (required)

### Profile configs
Nextflow provides common Nextflow pipeline configurations for comput environments [https://nf-co.re/configs/]

Some important aspects for configuration files include, resoruce allocaiton for your computing environment including but not limited to hpc configuration, container settings (conda, singlarity, and docker), and memory allocation.

Specifically for container configurations, I load miniconda for the conda package manager. Configuration settings for running -profile conda include:
```console
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
```
## MultiQC output 
In the docs folder [https://github.com/uel3/nf-UnO/tree/master/docs] there is a MultiQC Report for nf-uno analysis of the salmonella outbreak data set. 
Below is a breakdown of each section of the MultiQC Report:

### General Stats
General Stats is a summary table of several QC metrics for input reads before and after QC trimming and host removal. These metrics include % duplicates, total number of sequences (in millions), median seqeunce length, and % GC content.
![General Stats](https://github.com/uel3/nf-UnO/blob/master/images/general_stats.svg)

### MIDAS2 Species Abundance
MIDAS2 Species Abundance is a summary table of the taxonomic profile for input reads. These taxonomic profiles represent only those species with median coverage of at least 2 across the 15 SCGs and at least 50% of SCGs are uniquely mapped. The metrics in this table include:

Lineage: Taxonomic ID

Fraction Covered: Fraction of covered bases of reference in the MIDAS2 database (horizontal genome coverage).

Covered Bases: Number of bases covered by at least one post-filtered reads.

Genome Length: Length of the reference genome in MIDAS2 database.

Mean Coverage: Mean read depth across all covered bases (vertical genome coverage)

Total Depth: Total reads depth across all covered bases

Aligned Reads: Total read counts across covered bases before post-alignment filter

Mapped Reads: Total read count across covered bases after post-alignment filter (median marker coverage = 2, unique fraction covered= 0.5)

Continent: Source of reference genome in MIDAS2 database. 

![MIDAS2](https://github.com/uel3/nf-UnO/blob/master/images/NFUNOmidas2.svg)

### FastQC: Raw Reads
**spot for metrics to highlight/consider when evaluating how seqeuncing went**
![FastQC: Raw Reads](https://github.com/uel3/nf-UnO/blob/master/images/FastqcRaw.svg)

### Bowtie2: Host Removal 
Bowtie2 mapping statistics of input reads mapped to host for removal using GRCh38 reference genome. The user can view the mapping statistics as number of reads or percentage of total reads. The portion of mNGS reads in red are the only reads used in downstream processes, with human host sequences removed. mNGS reads that are blue, light blue, orange, or yellow are considered host reads are discarded, and not included in downstream analysis. 
![Host Removal](https://github.com/uel3/nf-UnO/blob/master/images/BT2HostRemoval.svg)

### FastQC: After Preprocessing
**spot for metrics to highlight/consider when evaluating how seqeuncing went but post preprocessing**
![FastQC: After Preprocessing](https://github.com/uel3/nf-UnO/blob/master/images/FastqcPostPreProcess.svg)

### Binning: Sample Coverage of DASTool Bins
Sample reads are mapped against DASTool refined bins to calculate bin depths with log transformation. 

Only bins with reads mapped from all samples are of interest for downstream analysis. nf-UnO can output DASTool refined bins only (seen in this example) or both refined DASTool bins and unrefined bins. If the user is interested in unrefined bins from the MaxBin2 and MetaBat2 binners along with DASTool refined bins, adding ```--postbinning_input 'both'``` to CLI pipeline execution will generate that output. 
From this heatmap output, the user can see the how the individual mNGS contribute to the bins generated by DASTool. The bins represent species identified from all outbreak mNGS reads. Bins with read depth values greater than -5.99 would indicate contribution by that individual mNGS read. Bins with read depth values equal to -6.0 would indicate no contribution by that individual mNGS read. Bins that have contributions from all individual mNGS reads will be considered potential etiological agents. Bins are oragnized in descending order by number of samples with non-zero coverage.
![Binning](https://github.com/uel3/nf-UnO/blob/master/images/Binning.svg)

### CheckM: Bin Quality Statistics
Qaulity assessement of genome bins using CheckM.

The user can see the quality metrics of the species bins generated from all outbreak reads. The user can interpret the quality of the species bins based on combinations of the following metrics: completeness, contamination, and strain heterogeneity. Genome size can also be indicative of species bin quality. Moderate to high quality bins will have >50% completeness and <10% contamination. Low quality bins will have completeness <50% and contamination <10%. Strain heterogeneity should be interpreted as the amount of contamination attributed to closely related species. Note the binning process is prone to binning closely related species together which can lead to contamination, but this type of contamination is not in any way related to mishandling in lab (ie contaminated specimens).
![CheckM](https://github.com/uel3/nf-UnO/blob/master/images/checkm.svg)

### nf-uno Workflow Summary
List of nextflow options including runName (nextflow assigned run name), launchDir (directory the pipeline was executed), workDir (directory where nf-uno processes are performed, these directories contain intermediate files that are not included in the output directory but may be of use to some users), projectDir (nf-UnO pipeline directory), userName, profile (package manager/container profiles used), configFiles (config file passed to the pipeline)

Input/output options including path to ```--input``` samplesheet.csv used to run nf-UnO and ```--outdir``` output directory for nf-UnO

Institutional config options including config_profile_description, config_profle_contact, and config_profile_url as part of institutional config used with ```-c```

### nf-uno Software Versions

### Software Versions
Versions of software used by nf-UnO
### nf-uno Methods Description
nf-UnO release version, Nextflow version, and execution command used to run nf-UnO. 

