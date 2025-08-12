process eggnog_mapper {
	container "quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"
	label "eggnog_mapper"
	label "process_high"
	tag "${speci}/${genome_id}"
	cpus 16
	time {4.d * task.attempt}
	memory {64.GB * task.attempt}

	input:
	tuple val(speci), val(genome_id), path(proteins)
	path(eggnog_db)

	output:
	tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.emapper.annotations"), emit: eggnog, optional: true
	tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.EMAPPER.DONE"), emit: done_sentinel

	script:
	"""
	set -e -o pipefail
	mkdir -p ${speci}/${genome_id}/	

	emapper.py -i ${proteins} --data_dir ${eggnog_db} --output ${speci}/${genome_id}/${genome_id} -m diamond --cpu $task.cpus --dmnd_algo 0
	touch ${speci}/${genome_id}/${genome_id}.emapper.annotations
	if [[ -z \$(grep -v "#" ${speci}/${genome_id}/${genome_id}.emapper.annotations) ]]; then
		rm ${speci}/${genome_id}/${genome_id}.emapper.annotations;
	fi
	touch ${speci}/${genome_id}/${genome_id}.EMAPPER.DONE
	"""

}
