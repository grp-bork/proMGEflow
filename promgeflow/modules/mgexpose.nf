// [speci, bin_id, gene_coords, conjscan, emapper, clusters, recombinases, genome_fa]
process mgexpose {
	publishDir path: params.output_dir, mode: "copy", pattern: "**.mge_islands.gff3", enabled: !params.tarball_output
	publishDir path: params.output_dir, mode: "copy", pattern: "**.mge_islands.ffn.gz", enabled: !params.tarball_output
	publishDir path: params.output_dir, mode: "copy", pattern: "**.gene_info.txt", enabled: !params.tarball_output
	label "annotate_genome"
	label "medium"
	label "mgexpose"
	container "ghcr.io/grp-bork/mgexpose:v3.10.0"
	// executor "local"  -> move to run.config @ EMBL
	tag "${speci}/${genome_id}"

	input:
	tuple val(speci), val(genome_id), path(gff), path(conjscan), path(emapper), path(gene_clusters), path(recombinases), path(genome_fa)
	path(mge_rules)
	path(conjscan_rules)
	path(phage_filter_terms)
	val(bulk_output)

	output:
	tuple val(speci), val(genome_id), path("**/*.mge_islands.gff3"), emit: gff, optional: true
	tuple val(speci), val(genome_id), path("**/*.unannotated_islands.gff3"), emit: genomic_islands, optional: true
	tuple val(speci), val(genome_id), path("**/*.NO_MGE"), emit: no_mge, optional: true
	tuple val(speci), val(genome_id), path("**/*.mge_islands.ffn.gz"), emit: fasta, optional: true
	tuple val(speci), val(genome_id), path("**/*.gene_info.txt"), emit: gene_info, optional: true
	tuple val(speci), val(genome_id), path("${genome_id}.pangenome.txt"), emit: pangenome_info, optional: true
	
	script:
	def y_cluster_option = (params.use_y_clusters) ? " --use_y_clusters --core_threshold -1" : ""
	def outdir = (bulk_output) ? "${speci}/${genome_id}" : "${genome_id}"

	"""
	mkdir -p ${outdir}/

	echo mgexpose denovo ${genome_id} \
			--input_genes ${gff} \
			--recombinases ${recombinases} \
			--speci ${speci} \
			--conjugation_data ${conjscan} \
			--phage_and_cargo_data ${emapper} \
			--cluster_data ${gene_clusters} \
			--output_dir ${outdir} \
			--genome_fasta ${genome_fa} \
			${y_cluster_option}
	mgexpose denovo ${genome_id} \
			--input_genes ${gff} \
			--recombinases ${recombinases} \
			--speci ${speci} \
			--conjugation_data ${conjscan} \
			--phage_and_cargo_data ${emapper} \
			--cluster_data ${gene_clusters} \
			--output_dir ${outdir} \
			--genome_fasta ${genome_fa} \
			${y_cluster_option}

	islands_gff=${outdir}/${genome_id}.mge_islands.gff3
	(grep mobile_genetic_element \${islands_gff} | grep -v mge= > ${genome_id}.NO_MGE) || true
	if [[ -s ${genome_id}.NO_MGE ]]; then 
		mv -v ${genome_id}.NO_MGE ${outdir}/
	else
		rm -f ${genome_id}.NO_MGE
	fi

	if [[ -f ${outdir}/${genome_id}.gene_info.txt ]]; then
		awk -v OFS='\\t' 'BEGIN {  n=core["True"]=core["False"]=0; } NR > 1 { core[\$10]++;n++;} END {print "#species","genome","n_genes","n_accessory","n_core","%acc"; printf("%s\\t%s\\t%s\\t%s\\t%s\\t%.2f\\n", "${speci}","${genome_id}",n,core["False"],core["True"],core["False"]/n*100);}' ${outdir}/${genome_id}.gene_info.txt > ${genome_id}.pangenome.txt
	fi

	rm -vf mgexpose.gff
	"""

}


process mgexpose_region {
	publishDir path: params.output_dir, mode: "copy", pattern: "**.mge_islands.gff3", enabled: !params.tarball_output
	publishDir path: params.output_dir, mode: "copy", pattern: "**.mge_islands.ffn.gz", enabled: !params.tarball_output
	publishDir path: params.output_dir, mode: "copy", pattern: "**.gene_info.txt", enabled: !params.tarball_output
	label "annotate_genome"
	label "mgexpose"
	container "ghcr.io/grp-bork/mgexpose:v3.10.0"
	// executor "local"  -> move to run.config @ EMBL
	tag "${speci}/${genome_id}"

	input:
	tuple val(speci), val(genome_id), val(region_id), path(gff), path(conjscan), path(emapper), path(recombinases), path(genome_fa)
	path(mge_rules)
	path(conjscan_rules)
	path(phage_filter_terms)
	val(bulk_output)

	output:
	tuple val(speci), val(genome_id), path("**/*.mge_islands.gff3"), emit: gff, optional: true
	tuple val(speci), val(genome_id), path("**/*.unannotated_islands.gff3"), emit: genomic_islands, optional: true
	tuple val(speci), val(genome_id), path("**/*.NO_MGE"), emit: no_mge, optional: true
	tuple val(speci), val(genome_id), path("**/*.mge_islands.ffn.gz"), emit: fasta, optional: true
	tuple val(speci), val(genome_id), path("**/*.gene_info.txt"), emit: gene_info, optional: true
	
	script:
	def outdir = (bulk_output) ? "${speci}/${genome_id}" : "${genome_id}"

	"""
	mkdir -p ${outdir}

	echo ${region_id} > region.txt

	echo mgexpose denovo ${genome_id} \
			--input_genes ${gff} \
			--recombinases ${recombinases} \
			--speci no_speci \
			--conjugation_data ${conjscan} \
			--phage_and_cargo_data ${emapper} \
			--contigs_are_islands \
			--output_dir ${outdir} \
			--genome_fasta ${genome_fa}

	mgexpose denovo ${genome_id} \
			--input_genes ${gff} \
			--recombinases ${recombinases} \
			--speci no_speci \
			--conjugation_data ${conjscan} \
			--phage_and_cargo_data ${emapper} \
			--contigs_are_islands \
			--output_dir ${outdir} \
			--genome_fasta ${genome_fa}

	islands_gff=${outdir}/${genome_id}.mge_islands.gff3
	(grep mobile_genetic_element \${islands_gff} | grep -v mge= > ${genome_id}.NO_MGE) || true
	if [[ -s ${genome_id}.NO_MGE ]]; then 
		mv -v ${genome_id}.NO_MGE ${outdir}/
	else
		rm -f ${genome_id}.NO_MGE
	fi
	"""

}
