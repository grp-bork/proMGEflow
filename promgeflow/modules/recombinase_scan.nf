process recombinase_scan {
	container "oras://ghcr.io/cschu/profile_me_ci:latest"
	tag "${genome_id}"
	cpus 4
	time {1.d * task.attempt}
	memory {8.GB * task.attempt}

	input:
	tuple val(speci), val(genome_id), path(proteins)
	path(recombinase_hmm_db)
	path(mge_rules)

	output:
	tuple val(speci), val(genome_id), path("${genome_id}/${genome_id}.recombinase_hmmsearch.out"), emit: recombinases_raw
	tuple val(speci), val(genome_id), path("${genome_id}/${genome_id}.recombinase_hmmsearch.besthits.out"), emit: recombinases, optional: true
	tuple val(speci), val(genome_id), path("${genome_id}/${genome_id}.recombinase_based_MGE_predictions.tsv"), emit: recomb_table, optional: true

	script:
	// gzip -dc ${genome_id}.faa.gz > ${genome_id}.faa
	"""
	mkdir -p ${genome_id}/
	
	hmmsearch -o /dev/null --cpu $task.cpus --tblout ${genome_id}/${genome_id}.recombinase_hmmsearch.out --cut_ga ${recombinase_hmm_db} ${proteins}
	parse_hmmsearch.py --mge_rules ${mge_rules} --prefix ${genome_id}/${genome_id} ${genome_id}/${genome_id}.recombinase_hmmsearch.out 
	if [[ ! -s ${genome_id}/${genome_id}.recombinase_hmmsearch.besthits.out ]]; then rm -rf ${genome_id}/${genome_id}.recombinase_hmmsearch.besthits.out; rm -rf ${genome_id}/${genome_id}.recombinase_based_MGE_predictions.tsv; fi
	"""
}
