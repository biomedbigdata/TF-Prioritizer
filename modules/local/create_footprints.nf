process CREATE_FOOTPRINTS {
    conda 'bioconda::bedtools'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3' :
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
    val searchType
    val maxDistance
    path peakFile

    output:
    path "footprints.bed" into footprints

    script:
    """
    
    """
}