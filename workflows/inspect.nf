/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowInspect.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { PEAKS } from '../subworkflows/local/peaks'
include { COUNTS } from '../subworkflows/local/counts'
include { DYNAMITE } from '../subworkflows/local/dynamite'
include { RANKING } from '../subworkflows/local/ranking'
include { CHIP_ATLAS } from '../subworkflows/local/chip_atlas'
include { EH_ATLAS } from '../subworkflows/local/eh_atlas'
include { FIMO } from '../subworkflows/local/fimo'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'
include { ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION as MAP_GTF } from '../modules/nf-core/atlasgeneannotationmanipulation/gtf2featureannotation/main.nf'
include { GAWK as REMOVE_CHR } from '../modules/nf-core/gawk/main'
include { GAWK as REMOVE_CHR_FAIDX } from '../modules/nf-core/gawk/main'
include { GAWK as CLEAN_BED } from '../modules/nf-core/gawk/main'
include { GAWK as CHR_M } from '../modules/nf-core/gawk/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { GAWK as CLEAN_INDEX } from '../modules/nf-core/gawk/main'
include { GAWK as REMOVE_GENE_VERSIONS } from '../modules/nf-core/gawk/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INSPECT {

    ch_versions = Channel.empty()

    ch_peaks = Channel.fromPath(params.input).splitCsv(header: true).map {
        entry ->
            checked_file = file(entry.file, checkIfExists: true);
            return [
                        [
                            id: entry.state + "_" + entry.antibody + "_" + checked_file.baseName,
                            state: entry.state,
                            antibody: entry.antibody
                        ],
                        checked_file
                    ]
    }
    ch_states = ch_peaks.map{meta, file -> meta.state}.toSortedList().flatten().unique()

    ch_state_pairs = ch_states  .combine(ch_states)
                                .filter{state1, state2 -> state1 < state2}

    MAP_GTF (
        [[id:"gtf"], params.gtf],
        [[], []]
    )

    REMOVE_GENE_VERSIONS(MAP_GTF.out.feature_annotation, [])

    ch_gene_map = REMOVE_GENE_VERSIONS.out.output

    SAMTOOLS_FAIDX (
        [[id:"genome"], params.fasta],
        [[], []]
    )

    REMOVE_CHR_FAIDX(SAMTOOLS_FAIDX.out.fai, [])
    // CHR_M(REMOVE_CHR_FAIDX.out.output, [])
    CLEAN_INDEX(REMOVE_CHR_FAIDX.out.output, [])

    CHIP_ATLAS(REMOVE_CHR_FAIDX.out.output, params.ehmm_antigens)
    EH_ATLAS(params.genome, params.tax_id, params.enhancer_atlas_tissues)
    CLEAN_BED([[id:"promoters"], params.promoters], [])
    REMOVE_CHR(CLEAN_BED.out.output, [])

    PEAKS(
        ch_peaks,
        ch_gene_map,
        CHIP_ATLAS.out.ehmm,
        EH_ATLAS.out.bed,
        REMOVE_CHR.out.output,
        CLEAN_INDEX.out.output,
    )

    ch_counts = Channel.value(file(params.rnaseq_counts, checkIfExists: true))
    ch_design = Channel.value(file(params.rnaseq_design, checkIfExists: true))

    COUNTS(ch_counts, ch_design, ch_gene_map)

    ch_affinityRatio_deseq = PEAKS.out.affinity_ratio.map{
        meta, file -> [meta.state1, meta.state2, meta, file]
    }.combine(COUNTS.out.deseq2.map{
        meta, file -> [meta.state1, meta.state2, meta, file]
    }, by: [0, 1]).map{
        state1, state2, meta1, affinity_ratio, meta2, deseq2 -> [meta1, affinity_ratio, deseq2]
    }

    DYNAMITE(
        ch_affinityRatio_deseq,
        params.dynamite_ofolds,
        params.dynamite_ifolds,
        params.dynamite_alpha,
        params.dynamite_randomize
    )

    ch_affinitySum_deseq_dynamite = PEAKS.out.affinity_sum.map{
        meta, file -> [meta.state1, meta.state2, meta.antibody, meta, file]
    }.combine(COUNTS.out.deseq2.map{
            meta, file -> [meta.state1, meta.state2, meta, file]
        }, by: [0, 1]).map{
            state1, state2, antibody, meta1, affinity_sum, meta2, deseq2 -> [state1, state2, antibody, meta1, affinity_sum, deseq2]
    }.combine(DYNAMITE.out.map{
            meta, file -> [meta.state1, meta.state2, meta.antibody, meta, file]
        }, by: [0, 1, 2]).map{
            state1, state2, antibody, meta1, affinity_sum, deseq2, meta2, dynamite -> [meta1, affinity_sum, deseq2, dynamite]
    }

    RANKING (
        ch_affinitySum_deseq_dynamite
    )

    FIMO(
        RANKING.out.tfs,
        PEAKS.out.stitched_enhancers,
        params.pwm
    )

    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
