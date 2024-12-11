process MIDAS2_HEATMAP {
    
    conda "conda-forge::python=3.9 conda-forge::pandas=1.3.0 anaconda::seaborn=0.11.0 conda-forge::pyyaml=6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:6687f364558e1e2f561f9867ee10f69e2cceff48-0' :
        'biocontainers/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:6687f364558e1e2f561f9867ee10f69e2cceff48-0' }"

    input:
    path combined_report
    tuple val(meta), path(reads)

    output:
    path "midas2_heatmap_mqc.txt", emit: midas2_mqc_heatmap
    path "versions.yml", emit: versions

    script:
    def samples = meta.collect{ it.id }.join(',')
    """
    midas2_heatmap_mqc.py -i $combined_report -s $samples -t 0.25
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        seaborn: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('seaborn').version)")
        pyyaml: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pyyaml').version)")
    END_VERSIONS
    """
}