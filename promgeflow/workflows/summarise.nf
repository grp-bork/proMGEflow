process pangenome_summary {
	executor "local"

	input:
	path(genome_report)
	path(speci_sizes)

	output:
	path("pangenome_summary.txt")

	script:
	"""
	head -n 1 ${genome_report} | awk -v OFS='\\t' '{ print \$0,"n_genomes"; }' >> pangenome_summary.txt
	join -1 1 -2 1 <(tail -n +2 ${genome_report} | sort -k1,1) ${speci_sizes} | sort -k1,1 -k2,2 | tr " " "\\t" >> pangenome_summary.txt
	"""
}


workflow summarise {
	take:
		genomes_ch
	main:
		genomes_ch.dump(pretty: true, tag: "summarise_input")

		genome_status_ch = genomes_ch.map { speci, genome_id, gdata, flags -> [ speci, genome_id, flags ] }
			
		genome_status_ch.dump(pretty: true, tag: "genome_status_ch")

		Channel.of(["#species", "genome", "has_genes", "has_species", "has_ref_clusters", "has_recombinases", "has_functional", "has_conjugation", "has_pangenome", "has_mges"])
			.concat(
				genome_status_ch
					.map { speci, genome_id, flags -> [
							speci, genome_id, flags.GENOME_ANNOTATION, flags.SPECIES_RECOGNITION, flags.SPECI_CLUSTER_SEQS, flags.RECOMBINASE_SCAN, flags.FUNCTIONAL_ANNOTATION, flags.CONJUGATION_SYSTEM_ANNOTATION, flags.PANGENOME_CLUSTERING, flags.MGE_ANNOTATION
						]
					}
			)
			.collectFile(name: "genome_status.txt", newLine: true, sort: true, storeDir: "${params.output_dir}") {
				item -> item.join("\t")
			}

		if (params.run_mode != "contig" && params.run_mode != "plasmid") {
			/* Generate a pangenome report for the input genomes with identifed specI */
			genome_summary_ch = genomes_ch
				.map { speci, genome_id, gdata, flags -> gdata.pangenome_info }
				.filter { it != null }
				.collectFile(name: "${workDir}/pangenome_info.txt", skip: 1, keepHeader: true, sort: true)

			pangenome_summary(genome_summary_ch, "${projectDir}/assets/speci_sizes_pg3.txt")
		}
}



	
