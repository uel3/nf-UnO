process EPI_REPORT {
    tag "${meta.id}"
    
    conda "conda-forge::pandas=2.0.0 conda-forge::numpy=1.24.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.0' :
        'quay.io/biocontainers/pandas:2.0.0' }"
        
    input:
    tuple val(meta), path(depth_data)
    path(checkm_data)
    path(gtdb_data)
    path(midas_data)
    
    output:
    tuple val(meta), path("*_mag_analysis.txt"), emit: report
    path "versions.yml"                        , emit: versions
    
    script:
    """
    epi_report.py \\
        --depth_data ${depth_data} \\
        --checkm_data ${checkm_data} \\
        --gtdb_data ${gtdb_data} \\
        --midas_data ${midas_data} \\
        --output ${meta.id}_epi_report.txt \\
        --html_output ${meta.id}_epi_report.html
        --min_completeness ${params.mag_analysis_min_completeness} \\
        --max_contamination ${params.mag_analysis_max_contamination} \\
        --min_ani ${params.mag_analysis_min_ani} \\
        --min_af ${params.mag_analysis_min_af}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}