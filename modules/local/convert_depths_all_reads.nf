process CONVERT_DEPTHS_ALL {
    tag "${meta.id}"
    conda "bioconda::bioawk=1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--hed695b0_5' :
        'biocontainers/bioawk:1.0--hed695b0_5' }"

    input:
    tuple val(meta), path(fasta), path(depth)
    //Adding empty val to maintain consitency with pipeline and force maxbin2 to read in the abund_list flag 
    output:
    tuple val(meta), path(fasta), val([]), val([]), path("*_abund_list.txt"), val([]), emit: output
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -f $depth

    # Determine the number of abundance columns
    n_abund=\$(awk 'NR==1 {print int((NF-3)/2)}' ${depth.toString() - '.gz'})

    # Generate abundance files for each read set
    for i in \$(seq 1 \$n_abund); do
        col=\$((i*2+2))
        bioawk -t '{if (NR > 1) {print \$1, \$'"\$col"'}}' ${depth.toString() - '.gz'} > ${prefix}_mb2_depth_\$i.txt
    done

    # Create a list of abundance files with full paths, each on a new line
    for file in ${prefix}_mb2_depth_*.txt; do
        echo "\$PWD/\$file" >> ${prefix}_abund_list.txt
    done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioawk: \$(bioawk --version | cut -f 3 -d ' ' )
    END_VERSIONS
    """
}