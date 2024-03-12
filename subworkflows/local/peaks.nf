include { BEDTOOLS_SORT as SORT_PEAKS } from "../../modules/nf-core/bedtools/sort/main"
include { BEDTOOLS_SORT as SORT_FOOTPRINTS } from "../../modules/nf-core/bedtools/sort/main"
include { MERGE_PEAKS } from "./merge_peaks"
include { BEDTOOLS_CLOSEST } from "../../modules/nf-core/bedtools/closest/main"
include { GAWK as FILTER_CLOSEST } from "../../modules/nf-core/gawk/main"
include { GAWK as CLEAN_BED_ANNOTATE_SAMPLES } from "../../modules/nf-core/gawk/main"
include { CAT_CAT as MERGE_ORIGINAL } from "../../modules/nf-core/cat/cat/main"
include { BEDTOOLS_MERGE } from "../../modules/nf-core/bedtools/merge/main"
include { BEDTOOLS_SUBTRACT as SUBTRACT_BLACKLIST } from "../../modules/nf-core/bedtools/subtract/main"
include { STARE } from "../../modules/local/stare/stare"
include { MERGE_IDENTICAL } from "../../modules/local/stare/merge_identical"
include { COMBINE_TABLES as AFFINITY_MEAN } from "../../modules/local/combine_tables"
include { COMBINE_TABLES as AFFINITY_SUM } from "../../modules/local/combine_tables"
include { COMBINE_TABLES as AFFINITY_RATIO } from "../../modules/local/combine_tables"
include { COPY } from "../../modules/local/copy"
include { EHMM } from './ehmm'
include { GAWK as REMOVE_CHR } from "../../modules/nf-core/gawk/main"
include { CHROMHMM } from "./chromhmm"
include { ROSE } from "./rose"

workflow PEAKS {
    take:
        ch_peaks
        annotation_map

        ch_background

        ch_enhancers
        ch_promoters

        ch_index

    main:
        CLEAN_BED_ANNOTATE_SAMPLES(ch_peaks, [])

        if (params.merge_peaks) {
            ch_sorted = MERGE_PEAKS(CLEAN_BED_ANNOTATE_SAMPLES.out.output).peaks
        } else {
            SORT_PEAKS (CLEAN_BED_ANNOTATE_SAMPLES.out.output, [])
            ch_sorted = SORT_PEAKS.out.sorted
        }

        if (params.peak_search_type == "inside") {
            ch_footprints = ch_sorted
        } else {
            COPY(ch_sorted.map{ meta, bed_file ->
                    return [meta, bed_file, bed_file.name + ".copy"]
                })

            BEDTOOLS_CLOSEST(COPY.out, [])

            FILTER_CLOSEST(BEDTOOLS_CLOSEST.out.output, [])

            if (params.peak_search_type == "incl_between") {
                ch_joined = FILTER_CLOSEST.out.output
                    .join(ch_sorted)
                    .map {
                        meta, closest, bed_file ->
                            [meta, [closest, bed_file]]
                    }

                MERGE_ORIGINAL(ch_joined)
                SORT_FOOTPRINTS(MERGE_ORIGINAL.out.file_out, [])
                BEDTOOLS_MERGE(SORT_FOOTPRINTS.out.sorted)
                ch_footprints = BEDTOOLS_MERGE.out.bed
            } else {
                ch_footprints = FILTER_CLOSEST.out.output
            }
        }

        REMOVE_CHR(ch_footprints, [])
        ch_footprints_nochr = REMOVE_CHR.out.output

        /*
        if (params.bam_design) {
            ch_bams = Channel.value(file(params.bam_design))
                                .splitCsv(header: true)
                                .map{ row -> [[ id: row["state"] + "_" + row["assay"], 
                                                state: row["state"], 
                                                assay: row["assay"]], 
                                                file(row["directory"])]}

            EHMM(
                ch_background,
                ch_enhancers,
                ch_promoters,
                ch_footprints_nochr,
                ch_bams,

                params.gtf,
                ch_index,

                params.ehmm_genomic_region_size,
                params.ehmm_train_split,
                params.ehmm_random_seed,
                params.ehmm_top_quantile,
                params.ehmm_n_samples,
                params.ehmm_n_states,
                params.ehmm_n_bins,
                params.ehmm_pseudocount
            )

            ch_footprints = EHMM.out.all
        }
        */

        if (params.bam_design2) {
            // Parse design file
            ch_bams = Channel.value(file(params.bam_design2))
                                .splitCsv(header: ["state", "assay", "bam", "control"], sep: '\t')
                                .map{
                                    row -> [
                                        [
                                            id: row["state"] + "_" + row["assay"],
                                            state: row["state"],
                                            assay: row["assay"],
                                        ],
                                        [
                                            file(row["bam"]), file(row["control"])
                                        ],
                                    ]
                                }
                                .multiMap{
                                    row -> 
                                        normal_bams:
                                            [row[0], row[1][0]]
                                        control_bams:
                                            [row[0], row[1][1]]
                                }
        // Reunite channels with unique controls
        ch_bams = ch_bams.normal_bams.mix(ch_bams.control_bams.unique{ it[1] })

        CHROMHMM(
            ch_bams,
            params.chromhmm_states,
            params.chromsizes
        )

        // Add id field to meta object for use in nf-core modules
        ch_chromhmm_enhancers = CHROMHMM.out
            .map{meta, path ->
                meta = meta + [id: "enhancers_" + meta.state]
                [meta, path]
            }
        
        ROSE(
            ch_chromhmm_enhancers,
            params.ucsc_file
        )

        // Adapt process output to previous design of ch_footprints
        ch_footprints = ROSE.out
            .map{meta, bed -> 
                meta = meta + [antibody: "enhancers"]
                [meta, bed]
            }

        }
        
        if (params.blacklist) {
            SUBTRACT_BLACKLIST(ch_footprints.map{ meta, bed_file -> [meta, bed_file, params.blacklist]})
            ch_blacklisted = SUBTRACT_BLACKLIST.out.bed
        } else {
            ch_blacklisted = ch_footprints
        }

        STARE (
            ch_footprints,
            params.pwm,
            params.gtf,
            params.fasta,
            params.stare_window_size,
            params.blacklist,
            params.stare_decay
        )

        ch_affinities = STARE.out.affinities

        MERGE_IDENTICAL(ch_affinities, annotation_map)

        if (params.merge_peaks)
        {
            ch_mean_affinities = MERGE_IDENTICAL.out
        } else {
            AFFINITY_MEAN(
                MERGE_IDENTICAL.out
                    .map{ meta, affinity_file ->
                            [
                                [
                                    id: meta.state + "_" + meta.antibody,
                                    state: meta.state,
                                    antibody: meta.antibody
                                ],
                                affinity_file
                            ]
                        }
                    .groupTuple(),
                "mean"
            )
            ch_mean_affinities = AFFINITY_MEAN.out
        }

        ch_affinity_antibody = ch_mean_affinities
            .map{ meta, affinity_file -> [meta, meta.antibody, affinity_file]}

        ch_affinity_pairings = ch_affinity_antibody
            .combine(
                ch_affinity_antibody,
                by: 1
            )
            .filter{
                antibody, meta1, affinites1, meta2, affinites2 ->
                    meta1.state < meta2.state
            }
            .map{
                antibody, meta1, affinites1, meta2, affinites2 ->
                    [[  id: meta1.state + ":" + meta2.state + "_" + antibody,
                        state1: meta1.state,
                        state2: meta2.state,
                        antibody: antibody
                    ],
                    [affinites1, affinites2]]
            }

        AFFINITY_RATIO(ch_affinity_pairings, "ratio")
        AFFINITY_SUM(ch_affinity_pairings, "sum")

    emit:
        affinity_ratio = AFFINITY_RATIO.out
        affinity_sum = AFFINITY_SUM.out
        peaks = ch_blacklisted
        stitched_enhancers = ROSE.out
}
