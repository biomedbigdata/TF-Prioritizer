process GET_OUTPUT {

	label "process_single"
    container "quay.io/biocontainers/pandas:1.4.3"

    input:
    tuple val(meta), path(emissions), path(bed)

    output:
	tuple val(meta), path("enhancers_${bed.baseName.split('_')[0]}.bed")

	script:
	"""
	get_chromhmm_results.py \
		--emissions $emissions \
		--bed $bed \
		--output enhancers_${bed.baseName.split('_')[0]}.bed
	"""

	stub:
	"""
	touch enhancers_${bed.baseName.split('_')[0]}.bed
	"""
}