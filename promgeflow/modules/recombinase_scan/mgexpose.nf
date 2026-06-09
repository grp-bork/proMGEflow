process mgexpose_recombinase_scan {
	container "ghcr.io/grp-bork/mgexpose:v3.10.0"
	tag "${speci}/${genome_id}"
	cpus 4
	time {1.d * task.attempt}
	memory {8.GB * task.attempt}
	label "small"
	label "mgexpose"

	input:
	tuple val(speci), val(genome_id), path(proteins), path(gff)
	path(recombinase_hmm_db)
	path(mge_rules)

	output:
	tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.recombinase_hmmsearch.out"), emit: recombinases_raw
	// tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.recombinase_hmmsearch.besthits.out"), emit: recombinases, optional: true
	// tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.recombinase_based_MGE_predictions.tsv"), emit: recomb_table, optional: true
	// tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.predicted_recombinase_mges.gff3"), emit: mge_pred_gff, optional: true
	tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.predicted_recombinase_mges.gff3"), emit: recombinases, optional: true
	tuple val(speci), val(genome_id), path("${speci}.${genome_id}.recombinase_sentinel"), emit: done_sentinel

	script:
	"""
	mkdir -p ${speci}/${genome_id}/

	mgexpose recombinase_scan -o ${speci}/${genome_id}/ -t ${task.cpus} ${proteins} ${gff} ${recombinase_hmm_db} ${genome_id}

	if [[ ! -s ${speci}/${genome_id}/${genome_id}.recombinase_hmmsearch.out ]]; then 
		rm -rvf ${speci}/${genome_id}/*		
	fi

	mv -v ${speci}/${genome_id}/${genome_id}.recombinase_scan.gff3 ${speci}/${genome_id}/${genome_id}.predicted_recombinase_mges.gff3 

	touch ${speci}.${genome_id}.recombinase_sentinel
	"""
	// parse_hmmsearch.py --proteins recombinase_scan.faa --mge_rules ${mge_rules} --prefix ${speci}/${genome_id}/${genome_id} ${speci}/${genome_id}/${genome_id}.recombinase_hmmsearch.out
}
