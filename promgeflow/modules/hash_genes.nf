process hash_genes {
	tag "${genome_id}"
	label "mgexpose"

	input:
	tuple val(genome_id), path(genes)

	output:
	tuple val(genome_id), path("${genome_id}.genes.sha256.tsv"), emit: genes_sha256

	script:
	"""
	mgexpose hash_genes -o . ${genes} ${genome_id} 
	"""

}