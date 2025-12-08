process macsyfinder {
	container "quay.io/biocontainers/macsyfinder:2.1.4--pyhdfd78af_0"
	tag "${speci}/${genome_id}"
	cpus 8
	time {1.d * task.attempt}
	memory {16.GB * task.attempt}
	label "small"

	input:
	tuple val(speci), val(genome_id), path(proteins)
	path(models)

	output:
	tuple val(speci), val(genome_id), path("**/${genome_id}.all_systems.tsv"), emit: macsy_report

	script:
	def prefix = (speci != null) ? "${speci}/${genome_id}" : "${genome_id}"

	"""
	set -e -o pipefail
	mkdir -p ${prefix}/

	if [[ "${proteins}" == *".gz" ]]; then
		gzip -dc ${proteins} > macsy.faa
	else
		ln -sf ${proteins} macsy.faa
	fi

	macsyfinder -vvv -w ${task.cpus} --models CONJ all --models-dir ${models} -o ${prefix} --db-type unordered --multi-loci all --sequence-db macsy.faa
	cp -v ${prefix}/all_systems.tsv ${prefix}/${genome_id}.all_systems.tsv || touch ${prefix}/${genome_id}.all_systems.tsv

	rm -vf macsy.faa
	"""
}
