## Parameters for current nf-UnO pipeline
### MIDAS2 default settings
MIDAS2 is set to run by default by `skip_midas2 = false`

midas2_run_snps uses the following threshold by default:

`midas2_snps_select_by = 'median_marker_coverage,unique_fraction_covered'`
`midas2_median_marker_coverage  = '2'`
`midas2_unique_fraction_covered = '0.5'`
### Trimmomatic default settings 

Trimmomatic adapter sequence is available in assets folder:

`adapter_seqeunce = "${projectDir}/assets/TruSeq3-PE.fa"`
`qual_trim = 20:30:10:8:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`
`method = PE`
### Host Removal default settings

Host removal default setting is set to `--host_genome 'GRCh38'`
### Metagenomic Co-assembly Settings/MEGAHIT default settings 

In nextflow.config `coassemble_group = true` is set to perform co-assembly on all samples with the same `group`, denoted in `--input` samplesheet.csv
In modules.config `ext.args = '--presets meta-large'` is set to account for complex metagenomes (i.e., stool specimens)

MEGAHIT memory settings are set to `--mem-flag 2` to use all memory to avoid timing out during coassembly step. The following settings are in place to allocate sufficient memory for MEGAHIT coassembly step but should be modified for different systems:

In the nextflow.config:

```console
megahit_fix_cpu_1          = false

max_memory                 = '128.GB'
max_cpus                   = 16

// Functions to fix number of cpus to allow reproducibility for MEGAHIT
// if corresponding parameters are specified, number of cpus is not increased with retries
def check_megahit_cpus (x, attempt ) {
    if (params.megahit_fix_cpu_1) return 1
    else return check_max (x * attempt, 'cpus' )
}
```
In base.config:
```console
withName: MEGAHIT {
        cpus          = { check_megahit_cpus (12, task.attempt  ) } 
        memory        = { params.max_memory }
        time          = { check_max (48.h  * task.attempt, 'time'   ) }
        errorStrategy = { task.exitStatus in [143,137,104,134,139,250] ? 'retry' : 'finish' }
    }
```

### Binning defaults settings 

Binning default settings are set to map reads to coassembly based on `group` information denoted in `--input` samplesheet.csv, post binning analyses will be performed on refined bins only.

In nextflow.config:
```console
binning_map_mode                     = 'group'
postbinning_input                    = 'refined_bins_only'
```
 