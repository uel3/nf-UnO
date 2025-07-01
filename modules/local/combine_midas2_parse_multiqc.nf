process COMBINE_MIDAS2_REPORTS {
    conda "conda-forge::python=3.9 conda-forge::pandas=1.1.5 conda-forge::pyyaml=6.*"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc4df9b7c2:ab110436faf952a33575c64dd74615a84011450b-0' :
        'quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc4df9b7c2:ab110436faf952a33575c64dd74615a84011450b-0' }"
    
    input:
    path snps_id_list

    output:
    path "combined_midas2_report_mqc.yaml", emit: combined_report
    path "versions.yml", emit: versions
    
    script:
    """
    combine_midas2_parse_mutliqc.py -i $snps_id_list -y combined_midas2_report_mqc.yaml
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        pyyaml: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pyyaml').version)")
    END_VERSIONS
    """
}

