include { GTDBTK_DB_PREPARATION } from '../../modules/local/gtdbtk_db_preparation'
include { GTDBTK_CLASSIFYWF     } from '../../modules/nf-core/gtdbtk/classifywf/main'

workflow GTDBTK {
    take:
    bins              // channel: [ val(meta), [bins] ]
    gtdb              // channel: path

    
    main:
    // Use collect() to get the file out of the channel
    if (gtdb.toString().endsWith('.tar.gz')) {
        // Expects to be tar.gz!
        ch_db_for_gtdbtk = GTDBTK_DB_PREPARATION(gtdb).db
    } else if (file(gtdb).isDirectory()) {
        // Directory handling - create a channel with the tuple
        gtdb_dir = file(gtdb).listFiles()
        ch_db_for_gtdbtk = Channel.value(["gtdb", gtdb_dir])
    } else {
        error("Unsupported object given to --gtdb, database must be supplied as either a directory or a .tar.gz file!")
    }
        
    GTDBTK_CLASSIFYWF (
        bins,
        ch_db_for_gtdbtk,
    )
    emit:
    summary     = GTDBTK_CLASSIFYWF.out.summary
    bac_summary = GTDBTK_CLASSIFYWF.out.bac_summary
    versions    = GTDBTK_CLASSIFYWF.out.versions
}