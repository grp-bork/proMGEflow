process linclust {
	container "quay.io/biocontainers/mmseqs2:13.45111--h2d02072_0"
	label "linclust"
	label "process_high"
	tag "${speci}/${genome_id}"
	memory {32.GB * task.attempt}
	time { 1.d * task.attempt }
	cpus 8

	input:
	tuple val(speci), val(genome_id), path(genes), path(speci_seqs)
	
	output:
	tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}_mmseqcluster.tsv.gz"), emit: mmseq_cluster
	tuple val(speci), val(genome_id), path("${speci}/${genome_id}/${genome_id}.LINCLUST_DONE"), emit: done_sentinel

	script:
	"""
	set -e -o pipefail
	mkdir -p tmp/ ${speci}/${genome_id}/

	gzip -c ${genes} > all_genes.fa.gz
	cat ${speci_seqs} >> all_genes.fa.gz
		
	mmseqs createdb all_genes.fa.gz mmseqsdb
	mmseqs linclust --threads ${task.cpus} --split-memory-limit ${task.memory.toGiga()}G --min-seq-id 0.95 -c 0.90 --cov-mode 0 mmseqsdb mmseqcluster ./tmp
	mmseqs createsubdb mmseqcluster mmseqsdb mmseqcluster_representatives
	mmseqs createtsv mmseqsdb mmseqsdb mmseqcluster mmseqcluster.tsv
	mv mmseqcluster.tsv ${speci}/${genome_id}/${genome_id}_mmseqcluster.tsv
	gzip -v ${speci}/${genome_id}/${genome_id}_mmseqcluster.tsv

	rm -rvf all_genes.fa.gz tmp/
	touch ${speci}/${genome_id}/${genome_id}.LINCLUST_DONE
	"""

}
