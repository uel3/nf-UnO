/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowUno.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.host_fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { INPUT_CHECK                } from '../subworkflows/local/input_check'
include { MIDAS2_DB                  } from '../subworkflows/local/midas2dbbuild'
include { MIDAS2_SPECIES_SNPS        } from '../modules/local/midas2/speciessnps'
include { BT2_HOST_REMOVAL_BUILD     } from '../modules/local/bowtie2/bt2_host_removal_build'
include { BT2_HOST_REMOVAL_ALIGN     } from '../modules/local/bowtie2/bt2_host_removal_align'
include { BINNING_PREP               } from '../subworkflows/local/binning_prep'
include { BINNING                    } from '../subworkflows/local/binning'
include { DASTOOL_BINNING_REFINEMENT } from '../subworkflows/local/dastool_binning_refinement'
include { DEPTHS                     } from '../subworkflows/local/depths'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_RAW                  } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED              } from '../modules/nf-core/fastqc/main'
include { TRIMMOMATIC                           } from '../modules/nf-core/trimmomatic/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { MEGAHIT                               } from '../modules/nf-core/megahit/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/* --  Create channel for host reference  -- */
if ( params.host_genome ) {
    host_fasta = params.genomes[params.host_genome].fasta ?: false
    ch_host_fasta = Channel
        .value(file( "${host_fasta}" ))
    host_bowtie2index = params.genomes[params.host_genome].bowtie2 ?: false
    ch_host_bowtie2index = Channel
        .value(file( "${host_bowtie2index}/*" ))
} else if ( params.host_fasta ) {
    ch_host_fasta = Channel
        .value(file( "${params.host_fasta}" ))
} else {
    ch_host_fasta = Channel.empty()
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow UNO {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )

    ch_raw_short_reads  = INPUT_CHECK.out.raw_short_reads
    ch_raw_long_reads   = INPUT_CHECK.out.raw_long_reads
    if ( !params.skip_midas2 ){
        MIDAS2_DB (
        
    )
        ch_versions = ch_versions.mix(MIDAS2_DB.out.midas2_db_version.first())
        MIDAS2_SPECIES_SNPS (MIDAS2_DB.out.midas2_db,
            ch_raw_short_reads
    )
    ch_versions = ch_versions.mix(MIDAS2_SPECIES_SNPS.out.versions.first())
    }
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //
    FASTQC_RAW (
        ch_raw_short_reads
    )
    TRIMMOMATIC {
        ch_raw_short_reads
    }
    ch_short_reads_prepped = Channel.empty()
    ch_short_reads_prepped = TRIMMOMATIC.out.trimmed_reads
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())
    //grouped_reads_ch = TRIMMOMATIC
        //.out
        //.trimmed_reads
        //.map { meta, reads ->
          //[ meta.group, meta, reads ]
        //}
        //.groupTuple(by: 0)
        //.map {
          //group, meta, reads ->

          //def groupedMeta = [:]
          //groupedMeta.id = "grouped_$group"
          //groupedMeta.group = group

          //def reads1 = reads.collect{ it[0] }
          //def reads2 = reads.collect{ it[1] }

          //[groupedMeta, reads1, reads2]
       // }
       //BT2_HOST_REMOVAL_BUILD only runs if a host_fasta is provided instead of host_genome
    if (params.host_fasta){
            BT2_HOST_REMOVAL_BUILD (
                ch_host_fasta
            )
            ch_host_bowtie2index = BT2_HOST_REMOVAL_BUILD.out.index
    }
    ch_bowtie2_removal_host_multiqc = Channel.empty()
    if (params.host_fasta || params.host_genome){
        BT2_HOST_REMOVAL_ALIGN (
            ch_short_reads_prepped,
            ch_host_bowtie2index
        )
        ch_short_reads_hostremoved = BT2_HOST_REMOVAL_ALIGN.out.reads
        ch_bowtie2_removal_host_multiqc = BT2_HOST_REMOVAL_ALIGN.out.log
        ch_versions = ch_versions.mix(BT2_HOST_REMOVAL_ALIGN.out.versions.first())
    } else {
        ch_short_reads_hostremoved = ch_short_reads_prepped
    }
    FASTQC_TRIMMED {
        ch_short_reads_hostremoved
    }
    ch_short_reads_assembly = Channel.empty()
    ch_short_reads_assembly = ch_short_reads_hostremoved
        .map {meta, reads ->
            def meta_new = meta - meta.subMap('run')
            [ meta_new, reads ]
        }
    
    /*
    CO-ASSEMBLY OF TRIMMED READS 
    */
    if (params.coassemble_group) {
            // short reads
            // group and set group as new id
        ch_short_reads_grouped = ch_short_reads_assembly
                    .map { meta, reads -> [ meta.group, meta, reads ] }
                    .groupTuple(by: 0)
                    .map { group, metas, reads ->
                        //params.bbnorm
                        def meta         = [:]
                        meta.id          = "group-$group"
                        meta.group       = group
                        [ meta, reads.collect { it[0] }, reads.collect { it[1] } ]
                    }
    } else {
            ch_short_reads_grouped = ch_short_reads_assembly
                .map { meta, reads -> [ meta, [ reads[0] ], [ reads[1] ] ] }
    }
            // long reads
            // group and set group as new id
    
    ch_assemblies = Channel.empty()
    MEGAHIT ( ch_short_reads_grouped )
            ch_megahit_assemblies = MEGAHIT.out.assembly
                .map { meta, assembly ->
                    def meta_new = meta + [assembler: 'MEGAHIT']
                    [ meta_new, assembly ]
                }
            ch_assemblies = ch_assemblies.mix(ch_megahit_assemblies)
            ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
    BINNING_PREP ( ch_assemblies, ch_short_reads_assembly )
            ch_versions = ch_versions.mix(BINNING_PREP.out.bowtie2_version.first())
    BINNING (BINNING_PREP.out.grouped_mappings, ch_short_reads_assembly)
        ch_bowtie2_assembly_multiqc = BINNING_PREP.out.bowtie2_assembly_multiqc
        ch_versions = ch_versions.mix(BINNING_PREP.out.bowtie2_version.first())
        ch_versions = ch_versions.mix(BINNING.out.versions)
    ch_binning_results_bins = BINNING.out.bins
    ch_binning_results_bins = ch_binning_results_bins
            .map { meta, bins ->
                def meta_new = meta + [refinement:'unrefined']
                [meta_new , bins]
            }
    ch_binning_results_unbins =  BINNING.out.unbinned
    ch_binning_results_unbins = ch_binning_results_unbins
            .map { meta, bins ->
                def meta_new = meta + [refinement:'unrefined_unbinned']
                [meta_new, bins]
            }
    ch_contigs_for_binrefinement = BINNING_PREP.out.grouped_mappings
                    .map{ meta, contigs, bam, bai -> [ meta, contigs ] }
    DASTOOL_BINNING_REFINEMENT ( ch_contigs_for_binrefinement, ch_binning_results_bins )
    ch_refined_bins = DASTOOL_BINNING_REFINEMENT.out.refined_bins
    ch_refined_unbins = DASTOOL_BINNING_REFINEMENT.out.refined_unbins
    ch_versions = ch_versions.mix(DASTOOL_BINNING_REFINEMENT.out.versions)
    //including the following channel mapping options in case we want to look at raw bins or both eventually
   if ( params.postbinning_input == 'raw_bins_only' ) {
        ch_input_for_postbinning_bins        = ch_binning_results_bins
        ch_input_for_postbinning_bins_unbins = ch_binning_results_bins.mix(ch_binning_results_unbins)
    } else if ( params.postbinning_input == 'refined_bins_only' ) {
        ch_input_for_postbinning_bins        = ch_refined_bins
        ch_input_for_postbinning_bins_unbins = ch_refined_bins.mix(ch_refined_unbins)
    } else if ( params.postbinning_input == 'both' ) {
        ch_all_bins = ch_binning_results_bins.mix(ch_refined_bins)
        ch_input_for_postbinning_bins        = ch_all_bins
        ch_input_for_postbinning_bins_unbins = ch_all_bins.mix(ch_binning_results_unbins).mix(ch_refined_unbins)
    } else {
        ch_input_for_postbinning_bins        = ch_binning_results_bins
        ch_input_for_postbinning_bins_unbins = ch_binning_results_bins.mix(ch_binning_results_unbins)
    }
    //map read depths to bins
    DEPTHS ( ch_input_for_postbinning_bins_unbins, BINNING.out.metabat2depths, ch_short_reads_assembly )
        ch_input_for_binsummary = DEPTHS.out.depths_summary
        ch_versions = ch_versions.mix(DEPTHS.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    
    // MODULE: MultiQC
    
    workflow_summary    = WorkflowUno.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowUno.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.raw_reads.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
