process JASPAR_DOWNLOAD {
    label 'process_single'

    container "registry.hub.docker.com/bigdatainbiomedicine/inspect-python"

    output:
        path("jaspar_motifs")

    script:
    """
        curl -o jaspar.zip https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_redundant_pfms_meme.zip
        unzip jaspar.zip -d jaspar_motifs
    """

    stub:
    """
    mkdir jaspar_motifs
    touch jaspar_motifs/MA0001.1.meme
    touch jaspar_motifs/MA0001.2.meme
    touch jaspar_motifs/MA0001.3.meme
    touch jaspar_motifs/MA0002.1.meme
    touch jaspar_motifs/MA0002.2.meme
    """
}