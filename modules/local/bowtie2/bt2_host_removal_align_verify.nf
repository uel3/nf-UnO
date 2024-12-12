/*
 * Bowtie2 for read removal
 */
process BT2_HOST_REMOVAL_ALIGN_VERIFY {
    tag "$meta.id"

    conda "bioconda::bowtie2=2.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.2--py38h1c8e9b9_1' :
        'biocontainers/bowtie2:2.4.2--py38h1c8e9b9_1' }"

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.verify.unmapped*.fastq.gz") , emit: reads_verify
    path  "*.verify.mapped*.read_ids.txt", optional:true , emit: read_ids_verify
    tuple val(meta), path("*.bowtie2_verify.log") , emit: log
    path "versions.yml"                           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def save_ids = (args2.contains('--host_removal_save_ids')) ? "Y" : "N"
    """
    bowtie2 -p ${task.cpus} \
            -x ${index[0].getSimpleName()} \
            -1 "${reads[0]}" -2 "${reads[1]}" \
            $args \
            --un-conc-gz ${prefix}.verify.unmapped_%.fastq.gz \
            --al-conc-gz ${prefix}.verify.mapped_%.fastq.gz \
            1> /dev/null \
            2> ${prefix}.bowtie2_verify.log
    if [ ! -f ${prefix}.verify.unmapped_1.fastq.gz ]; then
        cp "${reads[0]}" ${prefix}.verify.unmapped_1.fastq.gz
        cp "${reads[1]}" ${prefix}.verify.unmapped_2.fastq.gz
    fi
    if [ ${save_ids} = "Y" ] ; then
        gunzip -c ${prefix}.verify.mapped_1.fastq.gz | awk '{if(NR%4==1) print substr(\$0, 2)}' | LC_ALL=C sort > ${prefix}.verify.mapped_1.read_ids.txt
        gunzip -c ${prefix}.verify.mapped_2.fastq.gz | awk '{if(NR%4==1) print substr(\$0, 2)}' | LC_ALL=C sort > ${prefix}.verify.mapped_2.read_ids.txt
    fi
    rm -f ${prefix}.verify.mapped_*.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
    }