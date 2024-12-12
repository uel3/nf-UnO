include { BT2_HOST_REMOVAL_ALIGN        } from '../../modules/local/bowtie2/bt2_host_removal_align'
include { BT2_HOST_REMOVAL_ALIGN_VERIFY } from '../../modules/local/bowtie2/bt2_host_removal_align_verify'

workflow HOST_REMOVAL {
    take:
        ch_short_reads_prepped    // Input reads channel
        ch_host_bowtie2index    // Host genome index channel

    main:
        ch_versions = Channel.empty()
        // First host removal step
        BT2_HOST_REMOVAL_ALIGN(
            ch_short_reads_prepped,
            ch_host_bowtie2index
        )

        // Verification step
        BT2_HOST_REMOVAL_ALIGN_VERIFY(
            BT2_HOST_REMOVAL_ALIGN.out.reads,
            ch_host_bowtie2index
        )

    emit:
        reads = BT2_HOST_REMOVAL_ALIGN.out.reads         // Only emit reads from first removal
        removal_logs = BT2_HOST_REMOVAL_ALIGN.out.log    // Initial removal logs
        verify_logs = BT2_HOST_REMOVAL_ALIGN_VERIFY.out.log    // Verification logs
        versions = BT2_HOST_REMOVAL_ALIGN.out.versions.mix(BT2_HOST_REMOVAL_ALIGN_VERIFY.out.versions)
}
