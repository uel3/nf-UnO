//
// Check input samplesheet and get read channels
//

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    if(hasExtension(params.input, "csv")){
        try {
            // extracts read files from samplesheet CSV and distribute into channels
            ch_input_rows = Channel
                .from(file(params.input))
                .splitCsv(header: true)
                .map { row ->
                        if (row.size() < 4) {
                        error("Input samplesheet contains row with ${row.size()} column(s). Expects at least 4 (id, group, sr1, sr2).")
                        }
                            def id = row.sample ?: error("Sample ID is missing in the input samplesheet.")
                            //def run = row.run
                            def group = row.group ?: error("Group is missing in the input samplesheet.")
                            def sr1 = row.short_reads_1 ? file(row.short_reads_1, checkIfExists: true) : error("short_reads_1 is missing or invalid.")
                            def sr2 = row.short_reads_2 ? file(row.short_reads_2, checkIfExists: true) : error("sr2 is missing or invalid.")
                            def lr = row.long_reads ? file(row.long_reads, checkIfExists: true) : null
                            // Check if given combination is valid
                            //if (run != null && run == "") exit 1, "ERROR: Please check input samplesheet -> Column 'run' contains an empty field. Either remove column 'run' or fill each field with a value."
                            if (!sr1 || !sr2) {
                                error("Invalid input samplesheet: Both short_reads_1 and short_reads_2 must be provided for sample ${id}.")
                            }

                            [ id, group, sr1, sr2, lr ]
                }
        // Check for unique sample IDs
        ch_input_rows
            .map { it[0] }
            .toList()
            .map { ids ->
                if (ids.size() != ids.unique().size()) {
                    error("Duplicate sample IDs found in the input samplesheet. Each sample ID must be unique.")
                }
            }
        // separate short and long reads
        ch_raw_short_reads = ch_input_rows
            .map { id, group, sr1, sr2, lr ->
                        def meta = [:]
                        meta.id           = id
                        meta.group        = group
                            return [ meta, [ sr1, sr2 ] ]
                }
        ch_raw_long_reads = ch_input_rows
            .map { id, group, sr1, sr2, lr ->
                        if (lr) {
                            def meta = [:]
                            meta.id           = id
                            meta.group        = group
                            return [ meta, lr ]
                        } else {
                            return null
                        }
                    }
                    .filter { it != null }
        } catch (Exception e) {
            exit 1, "ERROR: Problem parsing input samplesheet: ${e.message}"
        }
    } else {
        ch_raw_short_reads = Channel
            .fromFilePairs(params.input, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
            .map { row ->
                        def meta = [:]
                        meta.id           = row[0]
                        meta.group        = 0
                        return [ meta, row[1] ]
            }
        ch_input_rows = Channel.empty()
        ch_raw_long_reads = Channel.empty()
    }

    // Ensure run IDs are unique within samples (also prevents duplicated sample names)

    // Note: do not need to check for PE/SE mismatch, as checks above do not allow mixing
    ch_input_rows
        .groupTuple(by: 0)
        .map { id, groups, sr1s, sr2s, lrs -> 
            if (groups.size() != groups.unique().size()) {
                error("ERROR: input samplesheet contains duplicated sample IDs! Check samplesheet for sample id: ${id}")
            }
        }
    emit:
    raw_short_reads  = ch_raw_short_reads
    raw_long_reads   = ch_raw_long_reads
    }