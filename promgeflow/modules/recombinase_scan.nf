process recombinase_scan {
	container "registry.git.embl.org/schudoma/recombinase_scan:v0.1"
	tag "${speci}/${genome_id}"
	cpus 4
	time {1.d * task.attempt}
	memory {8.GB * task.attempt}

	input:
	tuple val(speci), val(genome_id), path(proteins)
	path(recombinase_hmm_db)
	path(mge_rules)

	output:
	tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.recombinase_hmmsearch.out"), emit: recombinases_raw
	tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.recombinase_hmmsearch.besthits.out"), emit: recombinases, optional: true
	tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.recombinase_based_MGE_predictions.tsv"), emit: recomb_table, optional: true
	tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.predicted_recombinase_mges.gff3"), emit: mge_pred_gff, optional: true

	script:
	"""
	mkdir -p ${speci}/${genome_id}/
	
	hmmsearch -o /dev/null --cpu $task.cpus --tblout ${speci}/${genome_id}/${genome_id}.recombinase_hmmsearch.out --cut_ga ${recombinase_hmm_db} ${proteins}

	parse_hmmsearch.py --proteins ${proteins} --mge_rules ${mge_rules} --prefix ${speci}/${genome_id}/${genome_id} ${speci}/${genome_id}/${genome_id}.recombinase_hmmsearch.out

	if [[ ! -s ${speci}/${genome_id}/${genome_id}.recombinase_hmmsearch.besthits.out ]]; then 
		rm -rvf ${speci}/${genome_id}/${genome_id}.recombinase_hmmsearch.besthits.out;
		rm -rvf ${speci}/${genome_id}/${genome_id}.recombinase_based_MGE_predictions.tsv;
		rm -rvf ${speci}/${genome_id}/${genome_id}.predicted_recombinase_mges.gff3
	fi
	"""
}
