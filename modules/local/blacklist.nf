process BLACKLIST {
    conda 'bioconda::bedtools'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3' :
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
    path blacklist
    tuple val(sample), path(peakFile), val(hm), val(group)

    output:
    tuple val(sample), path("${peakFile.baseName}_blacklist-filtered.bed"), val(hm), val(group)

    script:
    """
    bedtools subtract -a "$peakFile" -b "$blacklist" -A > "${peakFile.baseName}_blacklist-filtered.bed" &&
    [ -s "${peakFile.baseName}_blacklist-filtered.bed" ]
    """
}