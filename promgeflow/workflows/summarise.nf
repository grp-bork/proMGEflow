include { publish_tarball; publish_gene_annotations; publish_recombinase_scan } from "../modules/publish"


process pangenome_summary {
	publishDir path: "${params.output_dir}", mode: "copy", enabled: !params.tarball_output
	label "tiny"
	label "summary"
	// executor "local"

	input:
	path(genome_report)
	path(speci_sizes)

	output:
	path("pangenome_summary.txt"), emit: "pangenome_summary"

	script:
	"""
	head -n 1 ${genome_report} | awk -v OFS='\\t' '{ print \$0,"n_genomes"; }' >> pangenome_summary.txt
	join -1 1 -2 1 <(tail -n +2 ${genome_report} | sort -k1,1) ${speci_sizes} | sort -k1,1 -k2,2 | tr " " "\\t" >> pangenome_summary.txt
	"""
}


process genome_status_summary {
	publishDir path: "${params.output_dir}", mode: "copy", enabled: !params.tarball_output
	label "tiny"
	label "summary"
	// executor "local"

	input:
	path(genome_data)

	output:
	path("genome_status_summary.txt"), emit: summary

	script:
	"""
	awk -v OFS='\\t' \\
		'BEGIN { print "#species", "genome", "has_genes", "has_species", "has_ref_clusters", "has_recombinases", "has_functional", "has_conjugation", "has_pangenome", "has_mges" } \\
		 { print \$0; }' ${genome_data} > genome_status_summary.txt
	"""
}


workflow summarise_and_publish {
	take:
		genomes_ch
	main:
		genomes_ch.dump(pretty: true, tag: "summarise_input")

		genome_status_ch = genomes_ch.map { speci, genome_id, gdata, flags -> [ speci, genome_id, flags ] }
			
		genome_status_ch.dump(pretty: true, tag: "genome_status_ch")

		results_ch = Channel.empty()

		if (params.run_mode != "contig" && params.run_mode != "plasmid") {
			/* Generate a pangenome report for the input genomes with identifed specI */
			genome_summary_ch = genomes_ch
				.map { speci, genome_id, gdata, flags -> gdata.pangenome_info }
				.filter { it != null }
				.collectFile(name: "${workDir}/pangenome_info.txt", skip: 1, keepHeader: true, sort: true)

			pangenome_summary(genome_summary_ch, "${projectDir}/assets/speci_sizes_pg3.txt")

			results_ch = results_ch.mix(pangenome_summary.out.pangenome_summary)
		}

		results_recombinases_ch = genomes_ch
			// .filter { it[2].recomb_table != null && it[2].recomb_gff != null && it[2].mge_gff == null }
			// .filter { it[2].recomb_gff != null && it[2].mge_gff == null }
			.filter { it[2].recombinases != null && it[2].mge_gff == null }
			// .map { speci, genome_id, gdata, flags -> [ speci, genome_id, [ gdata.proteins, gdata.genes, gdata.gff, gdata.recomb_table, gdata.recomb_gff ] ] }
			// .map { speci, genome_id, gdata, flags -> [ speci, genome_id, [ gdata.proteins, gdata.genes, gdata.gff, gdata.recomb_gff ] ] }
			.map { speci, genome_id, gdata, flags -> [ speci, genome_id, [ gdata.proteins, gdata.genes, gdata.gff, gdata.recombinases ] ] }
		
		results_mge_ch = genomes_ch
			.filter { it[2].mge_gff != null && it[2].mge_fasta != null }
			.map { speci, genome_id, gdata, flags -> [ speci, genome_id, [ gdata.proteins, gdata.genes, gdata.gff, gdata.mge_gff, gdata.mge_fasta ] ] }

		if (params.tarball_output) {
			results_ch = results_ch.mix(
				results_recombinases_ch
					.mix(results_mge_ch)
					.map { speci, genome_id, payload -> payload }
			)
			.collect()
			
			results_ch.dump(pretty: true, tag: "results_ch_sap")

			publish_tarball(results_ch, params.tarball_output)
			
		} else {

			genome_status_summary(genome_status_ch
				.map { speci, genome_id, flags -> ([speci, genome_id] + flags.values().collect { it -> String.valueOf(it) }).join("\t") }
				.collectFile(name: "genome_status_summary.txt.raw", newLine: true, sort: true, storeDir: "${workDir}" )
			)
			results_ch = results_ch.mix(genome_status_summary.out.summary)

			publish_gene_annotations(
				results_recombinases_ch.mix(results_mge_ch).map { speci, genome_id, payload -> [ speci, genome_id, [ payload[0], payload[1], payload[2] ] ] },
				params.simple_output
			)

			publish_recombinase_scan(
				results_recombinases_ch.map { speci, genome_id, payload -> [ speci, genome_id, payload[3], payload[4] ] },
				params.simple_output
			)

		}
}
