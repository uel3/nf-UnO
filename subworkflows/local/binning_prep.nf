// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { BOWTIE2_BUILDASSEMBLYINDEX } from '../../modules/local/bowtie2/buildassemblyindex'
include { BOWTIE2_ALIGNASSEMBLY      } from '../../modules/local/bowtie2/alignassembly'

workflow BINNING_PREP {

    take:
    assemblies           // channel: [ val(meta), path(assembly) ]
    reads                // channel: [ val(meta), [ reads ] ]

    main:

    BOWTIE2_BUILDASSEMBLYINDEX (assemblies)
    // build bowtie2 index from coassembly of all reads using combinations
    if(params.combinations){

        // read in combinations file
        combinations_file = file('/home/rsainsbury/workspace/uno/samplesheets/all_2-3set_combinations.txt')
        combinations_ch = Channel.fromPath(combinations_file)
            .splitText()
            .map { line ->
                def fields = line.trim().split(/\t/)
                def group = fields[0]
                def samples = fields[1].split(',') as List
                tuple(group, samples)
            }

        // remap the ch_short_reads_assembly channel so that that it is formatted correctly for joining to the combinations channel
        reads_by_sample_ch = reads.map { meta, reads ->
            tuple(meta["id"], meta, reads)
        }

        // flattening the combinations channel so that there is one row for each sample and group combination
        combinations_ch = combinations_ch.flatMap { group, samples ->
            samples.collect { sample ->
                tuple(sample, group)
            }
        }

        // combining the combinations channel with the reads by sample channel to get a channel where each row is a one group and pair of reads
        group_reads_ch = combinations_ch
            .combine(reads_by_sample_ch, by: 0)
            .map { sample, group, reads_meta, reads ->
                tuple(group, reads_meta, reads)
            }

        ch_bowtie2_input = BOWTIE2_BUILDASSEMBLYINDEX.out.bt2_index
            .map {meta, assembly, index -> [meta.group, meta, assembly, index ] }
            .combine(group_reads_ch, by:0)
            .map {group, assembly_meta, assembly, index, reads_meta, reads -> [assembly_meta, assembly, index, reads_meta, reads ]}

    // build bowtie2 index from coassembly of all reads using group 
    } else if (params.binning_map_mode == 'group'){
        // combine assemblies with reads of all samples
        ch_reads_bowtie2 = reads.map{ meta, reads -> [ meta.group, meta, reads] }
        ch_bowtie2_input = BOWTIE2_BUILDASSEMBLYINDEX.out.bt2_index
            .map {meta, assembly, index -> [meta.group, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by:0)
            .map {group, assembly_meta, assembly, index, reads_meta, reads -> [assembly_meta, assembly, index, reads_meta, reads ]}

    } else {
        // i.e. --binning_map_mode 'own'
        // combine assemblies (not co-assembled) with reads from own sample
        ch_reads_bowtie2 = reads.map{ meta, reads -> [ meta.id, meta, reads ] }
        ch_bowtie2_input = BOWTIE2_BUILDASSEMBLYINDEX.out.bt2_index
            .map { meta, assembly, index -> [ meta.id, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { id, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }
    }
    
    BOWTIE2_ALIGNASSEMBLY (ch_bowtie2_input)
    ch_grouped_mappings = BOWTIE2_ALIGNASSEMBLY.out.mappings
        .groupTuple(by:0)
        .map {meta, assembly, bams, bais -> [meta, assembly.sort()[0], bams, bais ] }

    emit:
    // TODO nf-core: edit emitted channels
    bowtie2_assembly_multiqc = BOWTIE2_ALIGNASSEMBLY.out.log.map { assembly_meta, reads_meta, log -> [ log ] }
    bowtie2_version          = BOWTIE2_ALIGNASSEMBLY.out.versions
    grouped_mappings         = ch_grouped_mappings
}