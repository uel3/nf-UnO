// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { MAXBIN2                              } from '../../modules/nf-core/maxbin2/main'
include { METABAT2_METABAT2                    } from '../../modules/nf-core/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { GUNZIP as GUNZIP_BINS                } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_UNBINS              } from '../../modules/nf-core/gunzip/main'

include { SPLIT_FASTA                          } from '../../modules/local/split_fasta'
//include { CONVERT_DEPTHS                       } from '../../modules/local/convert_depths'
include { CONVERT_DEPTHS_ALL                   } from '../../modules/local/convert_depths_all_reads'
include { ADJUST_MAXBIN2_EXT                   } from '../../modules/local/adjust_maxbin2_ext'

workflow BINNING {

    take:
    // TODO nf-core: edit input (take) channels
    assemblies           // channel: [ val(meta), path(assembly),path(bams), path(bais) ]
    reads                // channel: [ val(meta), [ reads ] ]
    //reads_list           // channel: [ val(meta), [ text_file_ch] ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    ch_summarizedepth_input = assemblies
                                .map { meta, assembly, bams, bais ->
                                        def meta_new = meta.clone()
                                    [ meta_new, bams, bais ]
                                }

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_summarizedepth_input )

    ch_metabat_depths = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .map { meta, depths ->
            def meta_new = meta.clone()
            meta_new['binner'] = 'MetaBAT2'

            [ meta_new, depths ]
        }

    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first())

    // combine depths back with assemblies
    ch_metabat2_input = assemblies
        .map { meta, assembly, bams, bais ->
            def meta_new = meta.clone()
            meta_new['binner'] = 'MetaBAT2'

            [ meta_new, assembly, bams, bais ]
        }
        .join( ch_metabat_depths, by: 0 )
        .map { meta, assembly, bams, bais, depths ->
            [ meta, assembly, depths ]
        }
    CONVERT_DEPTHS_ALL ( ch_metabat2_input )
        ch_maxbin2_input = CONVERT_DEPTHS_ALL.out.output
            .map { meta, assembly, reads_list, reads, abund_list, depth ->
                    def meta_new = meta + [binner: 'MaxBin2']
                [ meta_new, assembly, reads_list, reads, abund_list, depth ]
            }
        ch_versions = ch_versions.mix(CONVERT_DEPTHS_ALL.out.versions.first())
    // main bins for decompressing for MAG_DEPTHS
    ch_final_bins_for_gunzip = Channel.empty()
    // final gzipped bins
    ch_binning_results_gzipped_final = Channel.empty()
    // run binning
    METABAT2_METABAT2 ( ch_metabat2_input )
        // before decompressing first have to separate and re-group due to limitation of GUNZIP module
        ch_final_bins_for_gunzip = ch_final_bins_for_gunzip.mix( METABAT2_METABAT2.out.fasta.transpose() )
        ch_binning_results_gzipped_final = ch_binning_results_gzipped_final.mix( METABAT2_METABAT2.out.fasta )
        ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions.first())
    //adding a reads_list option to check binning perfomance without using precomputed depths
    //ch_maxbin2_input_readslist = assemblies
        //.map { meta, assembly, bams, bais ->
            //def meta_new = meta.clone()
            //meta_new['binner'] = 'MaxBin2'
            //[ meta_new, assembly ]
        //}
        //.combine(reads_list)
        //.map { meta, assembly, reads_list ->
            //[ meta, assembly, reads_list, [], [] ]
        //}
    MAXBIN2 ( ch_maxbin2_input )
    ADJUST_MAXBIN2_EXT ( MAXBIN2.out.binned_fastas )
        ch_final_bins_for_gunzip = ch_final_bins_for_gunzip.mix( ADJUST_MAXBIN2_EXT.out.renamed_bins.transpose() )
        ch_binning_results_gzipped_final = ch_binning_results_gzipped_final.mix( ADJUST_MAXBIN2_EXT.out.renamed_bins )
        ch_versions = ch_versions.mix(MAXBIN2.out.versions)
    ch_input_splitfasta = METABAT2_METABAT2.out.unbinned.mix(MAXBIN2.out.unbinned_fasta)
    SPLIT_FASTA ( ch_input_splitfasta )
    // large unbinned contigs from SPLIT_FASTA for decompressing for MAG_DEPTHS,
    // first have to separate and re-group due to limitation of GUNZIP module
    ch_split_fasta_results_transposed = SPLIT_FASTA.out.unbinned.transpose()
    ch_versions = ch_versions.mix(SPLIT_FASTA.out.versions)

    GUNZIP_BINS ( ch_final_bins_for_gunzip )
    ch_binning_results_gunzipped = GUNZIP_BINS.out.gunzip
        .groupTuple(by: 0)

    GUNZIP_UNBINS ( ch_split_fasta_results_transposed )
    ch_splitfasta_results_gunzipped = GUNZIP_UNBINS.out.gunzip
        .groupTuple(by: 0)

    ch_versions = ch_versions.mix(GUNZIP_BINS.out.versions.first())
    ch_versions = ch_versions.mix(GUNZIP_UNBINS.out.versions.first())
    emit:
    // TODO nf-core: edit emitted channels
    bins                                         = ch_binning_results_gunzipped
    bins_gz                                      = ch_binning_results_gzipped_final
    unbinned                                     = ch_splitfasta_results_gunzipped
    unbinned_gz                                  = SPLIT_FASTA.out.unbinned
    metabat2depths                               = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    versions                                     = ch_versions
}

