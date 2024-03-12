process FILTER_MOTIFS {
    container "registry.hub.docker.com/bigdatainbiomedicine/inspect-python"

    input:
        path(tfs_jaspar_ids)
        path(jaspar_motifs)
    
    output:
        path "sign_motifs/*.meme"
    
    script:
    """
        #!/usr/bin/env python3

        from os import mkdir
        from os.path import exists
        from shutil import copy

        tfs_jaspar_ids = "${tfs_jaspar_ids}"
        jaspar_motifs = "${jaspar_motifs}"

        # Read differentially expressed (DE) transcription factors (TF)
        with open(tfs_jaspar_ids, "r") as f:
            tfs_jaspar_ids = f.read().split('\\n')

        # Create directory for significant motif files
        mkdir("sign_motifs")

        # Iterate over TFs and store meme files for DE TFs
        for jaspar_id in tfs_jaspar_ids:
            if exists(f"jaspar_motifs/{jaspar_id}.meme"):
                copy(f"jaspar_motifs/{jaspar_id}.meme", f"sign_motifs/{jaspar_id}.meme")
    """

    stub:
    """
    touch motifs.meme
    """
}