process COLLECT_EXPRESSION {
    conda "pandas"
    container "tfprio-python"

    input:
        path(tpm)
        path(deseq)
        path(counts)
        val(ensgs)

    output:
        path("gene_expression")

    script:
        """
        collect_expression.py -t $tpm -c $counts -d $deseq -e ${ensgs.join(" ")} -o gene_expression
        """
}

process COLLECT_TF_DATA {
    container "tfprio-python"

    input:
        tuple val(symbol), path(biophysical_logo), path(model), path(jaspar_logos), val(ensg), path(heatmaps)

    output:
        tuple path("$symbol"), path("${symbol}.json")

    script:
        """
        mkdir -p $symbol
        ln -s ../$biophysical_logo $symbol/biophysical_logo.png
        ln -s ../$model $symbol/model.pwm
        ln -s ../$jaspar_logos $symbol/jaspar_logos
        ln -s ../$heatmaps $symbol/heatmaps

        create_directory_json.py -d $symbol -e $ensg -o ${symbol}.json
        """
}

process COLLECT_HEATMAPS {
    container 'ubuntu:22.04'

    input:
        path(heatmaps)

    output:
        path("*_heatmaps")

    script:
        """
        link_heatmaps.sh . .
        """
}

process COLLECT_RANKS {
    conda "pandas"
    container "tfprio-python"

    input:
        path(ranks)
    
    output:
        path("ranks.json")
    
    script:
        """
        collect_ranks.py -r ${ranks.join(' ')} -o ranks.json
        """
}

process COLLECT {
    conda "pandas"
    container "tfprio-python"

    input:
        path(gene_expression)
        path(ranks)
        path(tf_files)

    output:
        path("data")

    script:
        """
        mkdir -p data

        ln -s ../$ranks data/ranks.json
        ln -s ../$gene_expression data/gene_expression

        mkdir -p data/tfs
        
        for tf in $tf_files; do
            ln -s ../../\$tf data/tfs
        done
        """
}