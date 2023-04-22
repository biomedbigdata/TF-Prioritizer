#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnaseq
    Website: https://nf-co.re/rnaseq
    Slack  : https://nfcore.slack.com/channels/rnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta            = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.transcript_fasta = WorkflowMain.getGenomeAttribute(params, 'transcript_fasta')
params.additional_fasta = WorkflowMain.getGenomeAttribute(params, 'additional_fasta')
params.gtf              = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff              = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed         = WorkflowMain.getGenomeAttribute(params, 'bed12')
params.bbsplit_index    = WorkflowMain.getGenomeAttribute(params, 'bbsplit')
params.star_index       = WorkflowMain.getGenomeAttribute(params, 'star')
params.hisat2_index     = WorkflowMain.getGenomeAttribute(params, 'hisat2')
params.rsem_index       = WorkflowMain.getGenomeAttribute(params, 'rsem')
params.salmon_index     = WorkflowMain.getGenomeAttribute(params, 'salmon')
params.bwa_index        = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.bowtie2_index    = WorkflowMain.getGenomeAttribute(params, 'bowtie2')
params.chromap_index    = WorkflowMain.getGenomeAttribute(params, 'chromap')
params.blacklist        = WorkflowMain.getGenomeAttribute(params, 'blacklist')
params.macs_gsize       = WorkflowMain.getMacsGsize(params)
params.enable_conda     = false


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULE PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.hisat2_build_memory = params.rnaseq_hisat2_build_memory
params.fragment_size = params.chipseq_fragment_size
params.narrow_peak = params.chipseq_narrow_peak
params.min_reps_consensus = params.chipseq_min_reps_consensus
params.clip_r1 = params.chipseq_clip_r1
params.clip_r2 = params.chipseq_clip_r2
params.three_prime_clip_r1 = params.chipseq_three_prime_clip_r1
params.three_prime_clip_r2 = params.chipseq_three_prime_clip_r2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASEQ } from '../subworkflows/local/rnaseq'
include { PEAK_FILES } from '../subworkflows/local/peak_files'
include { INTEGRATE_DATA } from '../modules/local/dynamite'
include { PREPARE_FOR_CLASSIFICATION } from '../modules/local/dynamite'
include { DYNAMITE } from '../modules/local/dynamite'
include { DYNAMITE_FILTER } from '../modules/local/dynamite'
include { TFTG_SCORE } from '../modules/local/distributionAnalysis/tftg_scores'
include { RANKING } from '../modules/local/distributionAnalysis/ranking'
include { COLLECT_TFS } from '../modules/local/distributionAnalysis/collect_tfs'
include { TOP_TARGET_GENES } from '../modules/local/distributionAnalysis/top_target_genes'
include { BIOPHYSICAL_MODELS } from '../modules/local/biophysical_models'
include { HEATMAPS } from '../modules/local/distributionAnalysis/heatmaps'
include { TF_SEQUENCE } from '../modules/local/tf_sequence'

//
// WORKFLOW: Run main nf-core/rnaseq analysis pipeline
//
workflow TFPRIO {
    ch_versions = Channel.empty()

    RNASEQ ()
    PEAK_FILES ()

    ch_versions = ch_versions.mix(RNASEQ.out.versions, PEAK_FILES.out.versions)
    ch_affinity_ratios = PEAK_FILES.out.affinity_ratios
        .map { [it[1], it[2], it[0], it[3]] } // group1, group2, hm, affinityRatios
    ch_affinity_sums = PEAK_FILES.out.affinity_sums
        .map { [it[1], it[2], it[0], it[3]] } // group1, group2, hm, affinitySums
    ch_diff_expression = RNASEQ.out.deseq2 // group1, group2, deseq2-file
    ch_counts = RNASEQ.out.count // group1, group2, counts-file
    ch_rnaseq_samplesheet = RNASEQ.out.samplesheet // samplesheet-file

    ch_integration = ch_affinity_ratios
        .combine(ch_diff_expression, by: [0, 1]) // group1, group2, hm, affinityRatios, diffExpression
    
    INTEGRATE_DATA (ch_integration)
    PREPARE_FOR_CLASSIFICATION (INTEGRATE_DATA.out)
    DYNAMITE (
        PREPARE_FOR_CLASSIFICATION.out,
        Channel.value(params.dynamite_ofolds),
        Channel.value(params.dynamite_ifolds),
        Channel.value(params.dynamite_alpha),
        Channel.value(params.dynamite_randomize)
    )

    DYNAMITE_FILTER (
        DYNAMITE.out,
        Channel.value(params.dynamite_min_regression)
    )

    ch_tftg = DYNAMITE_FILTER.out
        .combine(ch_diff_expression, by: [0, 1]) // group1, group2, hm, coefficients, diffExpression
        .combine(ch_affinity_sums, by: [0, 1, 2]) // group1, group2, hm, coefficients, diffExpression, affinitySums
    
    TFTG_SCORE (ch_tftg) // group1, group2, hm, tftgScore
    RANKING (TFTG_SCORE.out) // group1, group2, hm, ranking

    COLLECT_TFS (
        RANKING.out.map { it[3] } .collect()
    ) // tfs

    TOP_TARGET_GENES (
        COLLECT_TFS.out,
        Channel.value(params.top_target_genes),
        ch_affinity_sums,
        ch_counts
    ) // group1, group2, hm, topTargetGenes

    ch_heatmaps = HEATMAPS (
        TOP_TARGET_GENES.out,
        ch_counts,
        ch_rnaseq_samplesheet,
        RNASEQ.out.ensg_map
    )   .transpose()
        .map { [it[0], it[1], it[2], it[3].name.replaceFirst(/^\w+:\w+_\w+_/, '').replaceFirst(/\.png$/, '') , it[3]] }
        // group1, group2, hm, tf, heatmap

    BIOPHYSICAL_MODELS (
        COLLECT_TFS.out,
        Channel.value(params.tepic_pwm)
    )

    ch_bio_logos = BIOPHYSICAL_MODELS.out.plots.flatten().map { [it.name.replaceAll(/.png$/, ''), it] }
    ch_bio_models = BIOPHYSICAL_MODELS.out.pwms.flatten().map { [it.name.replaceAll(/.pwm$/, ''), it] }

    ch_biophysical = ch_bio_logos
        .combine(ch_bio_models, by: 0) // tf, logo, model

    ch_jaspar = TF_SEQUENCE (COLLECT_TFS.out)
        .flatten()
        .map { [it.name.replaceAll(/_.*/, ''), it] }
        .groupTuple(by: 0) // tf, [logos]

    ch_logos = ch_biophysical
        .combine(ch_jaspar, by: 0) // tf, biophysical_logo, model, [jaspar_logos]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    TFPRIO ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/