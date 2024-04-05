process MIDAS2_SPECIES_SNPS {
    //errorStrategy 'ignore'
    tag "$meta.id"
    label 'process_medium'

    conda {
        name = 'midas2'
        file = '/scicomp/home-pure/uel3/UnO_nf/nf-core-uno/modules/local/midas2/environment.yml'
    }

    input:
    tuple val(meta), path(reads)

    output:
    path( "midas2_output/${meta.id}/species/log.txt" )
    path( "midas2_output/${meta.id}/species/species_profile.tsv" ), emit: species_id
    path( "midas2_output/${meta.id}/temp/*" ), optional: true //adding the optional: true keeps nf from throwing error
    path( "midas2_output/${meta.id}/snps/log.txt" )
    path( "midas2_output/${meta.id}/snps/snps_summary.tsv"), emit: midas2_snps
    path( "midas2_output/${meta.id}/snps/*.snps.tsv.lz4" )
    path( "midas2_output/${meta.id}/bt2_indexes/snps/*" ), optional: true
   

    script: //getting an error that midas2 cannot find hs-blastn but it is in the midas_changes env located :/scicomp/home-pure/uel3/.conda/envs/midas_changed/bin/hs-blastn
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def midas2_dbname = "--midasdb_name uhgg"
    def midas2_dbdir = "--midasdb_dir my_midasdb_uhgg"
    def outdir = "midas2_output"
    """
    midas2 run_species \
      --sample_name $prefix \
      -1 ${reads[0]} \
      -2 ${reads[1]} \
      $midas2_dbname \
      $midas2_dbdir \
      --num_cores $task.cpus \
      $outdir
    midas2 run_snps \
      --sample_name $prefix \
      -1 ${reads[0]} \
      -2 ${reads[1]} \
      $midas2_dbname \
      $midas2_dbdir \
      $args \
      --num_cores $task.cpus \
      $outdir
    """

    stub:
    """
    mkdir midas2_output
    mkdir midas2_output/${meta.id}
    mkdir midas2_output/${meta.id}/species
    touch midas2_output/${meta.id}/species/log.txt
    touch midas2_output/${meta.id}/species/species_profile.tsv
    mkdir midas2_output/${meta.id}/temp
    touch midas2_output/${meta.id}/temp/stub
    mkdir midas2_output/${meta.id}/bt2_indexes
    mkdir midas2_output/${meta.id}/bt2_indexes/snps
    touch midas2_output/${meta.id}/bt2_indexes/snps/stub
    mkdir midas2_output/${meta.id}/snps
    touch midas2_output/${meta.id}/snps/log.txt
    touch midas2_output/${meta.id}/snps/snps_summary.tsv
    touch midas2_output/${meta.id}/snps/stub.snps.tsv.lz4
    """
}