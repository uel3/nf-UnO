include { GTDBTK_DB_PREPARATION } from '../../modules/local/gtdbtk_db_preparation'
include { GTDBTK_CLASSIFYWF     } from '../../modules/nf-core/gtdbtk/classifywf/main'

workflow GTDBTK {
    take:
    bins              // channel: [ val(meta), [bins] ]
    gtdb              // channel: path
    
    main:
    if ( gtdb.extension == 'gz' ) {
        // Expects to be tar.gz!
            ch_db_for_gtdbtk = GTDBTK_DB_PREPARATION ( gtdb ).db
    } else if ( gtdb.isDirectory() ) {
        // The classifywf module expects a list of the _contents_ of the GTDB
        // database, not just the directory itself (I'm not sure why). But
        // for now we generate this list before putting into a channel,
        // then grouping again to pass to the module.
        // Then make up meta id to match expected channel cardinality for GTDBTK
        gtdb_dir = gtdb.listFiles()
        ch_db_for_gtdbtk = Channel
                            .of(gtdb_dir)
                            .collect()
                            .map { ["gtdb", it] }
    } else {
        error("Unsupported object given to --gtdb, database must be supplied as either a directory or a .tar.gz file!")
    }
    GTDBTK_CLASSIFYWF (
        ch_filtered_bins.passed.groupTuple(),
        ch_db_for_gtdbtk,
        params.gtdbtk_pplacer_useram ? false : true,
        gtdb_mash
    )
    emit:
    summary     = GTDBTK_CLASSIFYWF.out.summary
    versions    = GTDBTK_CLASSIFYWF.out.versions
}