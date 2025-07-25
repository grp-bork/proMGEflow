process eggnog_mapper {
	container "quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"
	label "eggnog_mapper"
	tag "${genome_id}"
	cpus 16
	time {4.d * task.attempt}
	memory {64.GB * task.attempt}

	input:
	tuple val(speci), val(genome_id), path(proteins)
	path(eggnog_db)

	output:
	tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.emapper.annotations"), emit: eggnog

	script:
	"""
	mkdir -p ${speci}/${genome_id}/	

	emapper.py -i ${proteins} --data_dir ${eggnog_db} --output ${speci}/${genome_id}/${genome_id} -m diamond --cpu $task.cpus --dmnd_algo 0
	"""

}
