process RUN_FIMO {
    tag "${meta.motif}"

    container "registry.hub.docker.com/leonhafner/memesuite"

    input:
        tuple val(meta), path(motif_file)
        path sequence_file
    
    output:
        tuple val(meta), path("fimo_${meta.motif}")
    
    script:
    """
    fimo --o fimo_${meta.motif} --max-stored-scores 1000000 ${motif_file} ${sequence_file}
    """

    stub:
    """
    mkdir fimo_${meta.motif}
    touch fimo_${meta.motif}/best_site.narrowPeak
    touch fimo_${meta.motif}/cisml.xml
    touch fimo_${meta.motif}/fimo.gff
    touch fimo_${meta.motif}/fimo.html
    touch fimo_${meta.motif}/fimo.tsv
    touch fimo_${meta.motif}/fimo.xml
    """

    
}