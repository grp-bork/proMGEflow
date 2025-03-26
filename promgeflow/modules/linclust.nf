process linclust {
	container "quay.io/biocontainers/mmseqs2:13.45111--h2d02072_0"
	label "linclust"
	tag "${genome_id}"
	memory {32.GB * task.attempt}
	time { 1.d * task.attempt }
	cpus 8

	input:
	tuple val(speci), val(genome_id), path(genes), path(speci_seqs)
	
	output:
	tuple val(speci), val(genome_id), path("${genome_id}/${genome_id}_mmseqcluster.tsv.gz"), emit: mmseq_cluster
	// tuple val(speci), val(genome_id), path("${genome_id}/${genome_id}_speci_headers.gz"), emit: speci_headers
	// tuple val(speci), val(genome_id), path("${genome_id}/${genome_id}_gene_headers.gz"), emit: gene_headers
	tuple val(speci), val(genome_id), path("${genome_id}/${genome_id}.DONE"), emit: done_sentinel

	script:
	// gzip -dc ${speci_seqs} | sed "s/^>/>proMGE__/" | gzip -c - > all_genes.fa.gz  // potentially unnecessary
	"""
	set -e -o pipefail
	mkdir -p tmp/ ${genome_id}/

	gzip -c ${genes} > all_genes.fa.gz
	cat ${speci_seqs} >> all_genes.fa.gz
		
	mmseqs createdb all_genes.fa.gz mmseqsdb
	mmseqs linclust --threads ${task.cpus} --split-memory-limit 150G --min-seq-id 0.95 -c 0.90 --cov-mode 0 mmseqsdb mmseqcluster ./tmp
	mmseqs createsubdb mmseqcluster mmseqsdb mmseqcluster_representatives
	mmseqs createtsv mmseqsdb mmseqsdb mmseqcluster mmseqcluster.tsv
	mv mmseqcluster.tsv ${genome_id}/${genome_id}_mmseqcluster.tsv
	gzip ${genome_id}/${genome_id}_mmseqcluster.tsv

	rm -rvf all_genes.fa.gz tmp/
	touch ${genome_id}/${genome_id}.DONE
	"""
	// gzip -dc ${genes} | grep "^>" | cut -f 1 -d " " | sed "s/^>//" | gzip -c - > ${genome_id}/${genome_id}_gene_headers.gz

	// gzip -dc all_genes.fa.gz | grep "^>" | cut -f 1 -d " " | sed "s/^>//" | gzip -c - > ${genome_id}/${genome_id}_cluster_input.gz

	// gzip -dc ${speci_seqs} | grep "^>" | cut -f 1 -d " " | sed "s/^>//" | gzip -c - > ${genome_id}/${genome_id}_speci_headers.gz
}
