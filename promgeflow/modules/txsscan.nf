process txsscan {
	container = "quay.io/biocontainers/macsyfinder:2.1.4--pyhdfd78af_0"
	tag "${genome_id}"
	cpus 8
	time {1.d * task.attempt}
	memory {16.GB * task.attempt}

	input:
	tuple val(speci), val(genome_id), path(proteins)
	path(txsscan_models)

	output:
	tuple val(speci), val(genome_id), path("**/all_systems.tsv"), emit: txsscan_report

	script:
	def prefix = (speci != null) ? "${speci}/${genome_id}" : "${genome_id}"

	"""
	mkdir -p ${prefix}/
	macsyfinder -vvv -w ${task.cpus} --models TXSS all --models-dir ${txsscan_models} -o ${prefix} --db-type unordered --multi-loci all --sequence-db ${proteins}
	"""
}
