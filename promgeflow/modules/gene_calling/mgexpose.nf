process mgexpose_call_genes {
	tag (batch_size > 1) ? "batch_${batch_id}" : genomes[0]
	label "small"
	label "mgexpose"
	container "ghcr.io/grp-bork/mgexpose:v3.10.0"
	cpus 1
	time {1.d * task.attempt}
	memory {8.GB * task.attempt}

	input:
	tuple val(batch_id), path(genomes)
	val(batch_size)

	output:
	path("gene_calls/**.*"), emit: annotations

	script:
	"""
	mkdir -p gene_calls/
	for genome_file in ${genomes}; do
		genome_id=\$(basename \$genome_file)
		echo \$genome_id
		
		# prodigal -i \$(basename \$genome_file .gz) -f gff -o prodigal/\$genome_id.gff -a prodigal/\$genome_id.faa -d prodigal/\$genome_id.ffn
		mgexpose call_genes -o gene_calls/\$genome_id -t ${task.cpus} \$genome_file \$genome_id
		
	done
	"""

}
