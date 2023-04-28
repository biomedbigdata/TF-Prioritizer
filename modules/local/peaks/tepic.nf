process TEPIC {
    label 'process_medium'

    conda "bioconda::bedtools conda-forge::gxx conda-forge::r-base conda-forge::r-gplots conda-forge::r-ggplot2 conda-forge::r-glmnet conda-forge::r-doMC conda-forge::r-reshape2 conda-forge::r-gridExtra pandas"
    container "tepic"

    input:
    tuple val(hm), val(group), val(sample), path(experimental_bed)
    path(pwm)
    path(gtf)
    path(reference_fasta)
    val(windowSize)
    val(loopWindows)
    val(exponentialDecay)
    val(normalizePeakLength)
    val(maxMinutesPerChromosome)
    val(originalScaling)
    val(pValue)

    output:
    tuple val(hm), val(group), path("*_sequences.tsv"), emit: sequences
    tuple val(hm), val(group), path("*_Affinity_Gene_View_Filtered.txt"), emit: affinities
    tuple val(sample), val(hm), val(group), path("*_canidate_binding_regions_Filtered_Regions.bed"), emit: filtered_regions

    // TODO: Implement missing config parameters

    script:
    """
    g++ $projectDir/assets/tepic/code/TRAPmulti.cpp -O3 -fopenmp -o $projectDir/assets/tepic/code/TRAPmulti

    bash $projectDir/assets/tepic/code/TEPIC.sh \
        -b $experimental_bed \
        -o $sample \
        -S ${sample}_sequences.tsv \
        -g $reference_fasta \
        -a $gtf \
        -f $gtf \
        -p $pwm \
        -w $windowSize \
        -s $loopWindows \
        -l ${!normalizePeakLength} \
        -i $maxMinutesPerChromosome \
        -v $pValue \
        -c $task.cpus
    """
}