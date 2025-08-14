// [speci, bin_id, gene_coords, txsscan, emapper, clusters, recombinases, genome_fa]
process mgexpose {
	label "annotate_genome"
	label "medium"
	container "ghcr.io/cschu/mgexpose:v3.7.6"
	executor "local"
	tag "${speci}/${genome_id}"

	input:
	tuple val(speci), val(genome_id), path(gff), path(txsscan), path(emapper), path(gene_clusters), path(recombinases), path(genome_fa)
	path(mge_rules)
	path(txsscan_rules)
	path(phage_filter_terms)
	val(simple_output)

	output:
	tuple val(speci), val(genome_id), path("**/*.mge_islands.gff3"), emit: gff, optional: true
	path("**/*.unannotated_islands.gff3"), emit: genomic_islands, optional: true
	path("**/*.NO_MGE"), emit: no_mge, optional: true
	path("**/*.mge_islands.ffn.gz"), emit: fasta, optional: true
	path("**/*.gene_info.txt"), emit: gene_info, optional: true
	
	script:
	def y_cluster_option = (params.use_y_clusters) ? " --use_y_clusters" : ""
	def outdir = (simple_output) ? "${genome_id}" : "${speci}/${genome_id}"

	"""
	mkdir -p ${outdir}/
	echo mgexpose denovo ${genome_id} ${gff} ${recombinases} ${mge_rules} \
			--speci ${speci} \
			--txs_macsy_rules ${txsscan_rules} \
			--txs_macsy_report ${txsscan} \
			--phage_eggnog_data ${emapper} \
			--phage_filter_terms ${phage_filter_terms} \
			--cluster_data ${gene_clusters} \
			--output_dir ${outdir} \
			--write_gff \
			--write_genes_to_gff \
			--add_functional_annotation \
			--dump_genomic_islands \
			--extract_islands ${genome_fa} \
			--output_suffix mge_islands \
			${y_cluster_option}
	mgexpose denovo ${genome_id} ${gff} ${recombinases} ${mge_rules} \
			--speci ${speci} \
			--txs_macsy_rules ${txsscan_rules} \
			--txs_macsy_report ${txsscan} \
			--phage_eggnog_data ${emapper} \
			--phage_filter_terms ${phage_filter_terms} \
			--cluster_data ${gene_clusters} \
			--output_dir ${outdir} \
			--write_gff \
			--write_genes_to_gff \
			--add_functional_annotation \
			--dump_genomic_islands \
			--extract_islands ${genome_fa} \
			--output_suffix mge_islands \
			${y_cluster_option}

	islands_gff=${outdir}/${genome_id}.mge_islands.gff3
	(grep mobile_genetic_element \${islands_gff} | grep -v mge= > ${genome_id}.NO_MGE) || true
	if [[ -s ${genome_id}.NO_MGE ]]; then 
		mv -v ${genome_id}.NO_MGE ${outdir}/
	else
		rm -f ${genome_id}.NO_MGE
	fi
	"""

}


process mgexpose_region {
	label "annotate_genome"
	container "ghcr.io/cschu/mgexpose:v3.7.6"
	executor "local"
	tag "${speci}/${genome_id}"

	input:
	tuple val(speci), val(genome_id), path(gff), path(txsscan), path(emapper), val(region_id), path(recombinases), path(genome_fa)
	path(mge_rules)
	path(txsscan_rules)
	path(phage_filter_terms)

	output:
	tuple val(speci), val(genome_id), path("**/*.mge_islands.gff3"), emit: gff, optional: true
	path("**/*.unannotated_islands.gff3"), emit: genomic_islands, optional: true
	path("**/*.NO_MGE"), emit: no_mge, optional: true
	path("**/*.mge_islands.ffn.gz"), emit: fasta, optional: true
	path("**/*.gene_info.txt"), emit: gene_info, optional: true
	
	script:
	def y_cluster_option = (params.use_y_clusters) ? " --use_y_clusters" : ""
	def outdir = "${speci}/${genome_id}"

	"""
	mkdir -p ${outdir}

	echo ${region_id} > region.txt

	echo mgexpose denovo ${genome_id} ${gff} ${recombinases} ${mge_rules} \
			--speci no_speci \
			--txs_macsy_rules ${txsscan_rules} \
			--txs_macsy_report ${txsscan} \
			--phage_eggnog_data ${emapper} \
			--phage_filter_terms ${phage_filter_terms} \
			--precomputed_islands region.txt \
			--output_dir ${outdir} \
			--write_gff \
			--write_genes_to_gff \
			--add_functional_annotation \
			--dump_genomic_islands \
			--extract_islands ${genome_fa} \
			--output_suffix mge_islands

	mgexpose denovo ${genome_id} ${gff} ${recombinases} ${mge_rules} \
			--speci no_speci \
			--txs_macsy_rules ${txsscan_rules} \
			--txs_macsy_report ${txsscan} \
			--phage_eggnog_data ${emapper} \
			--phage_filter_terms ${phage_filter_terms} \
			--precomputed_islands region.txt \
			--output_dir ${outdir} \
			--write_gff \
			--write_genes_to_gff \
			--add_functional_annotation \
			--dump_genomic_islands \
			--extract_islands ${genome_fa} \
			--output_suffix mge_islands

	islands_gff=${outdir}/${genome_id}.mge_islands.gff3
	(grep mobile_genetic_element \${islands_gff} | grep -v mge= > ${genome_id}.NO_MGE) || true
	if [[ -s ${genome_id}.NO_MGE ]]; then 
		mv -v ${genome_id}.NO_MGE ${outdir}/
	else
		rm -f ${genome_id}.NO_MGE
	fi
	"""

}
