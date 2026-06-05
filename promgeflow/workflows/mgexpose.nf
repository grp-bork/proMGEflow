include { mgexpose; mgexpose_region } from "../modules/mgexpose"
include { publish_gene_annotations } from "../modules/publish"

workflow mgexpose_denovo {
	take:
		genomes_ch
	main:

		mge_ch = genomes_ch()
		if (params.run_mode == "contig" || params.run_mode == "plasmid") {
			annotation_data_ch = genomes_ch
				.map { speci, genome_id, gdata -> [ speci, genome_id, "dummy_region", gdata.gff, gdata.conjugation_system_data, gdata.emapper, gdata.recombinases, gdata.genome ] }
			annotation_data_ch.dump(pretty: true, tag: "annotation_data_ch")

			mgexpose_region(
				annotation_data_ch,
				"${projectDir}/assets/mge_rules_ms.txt",
				"${projectDir}/assets/conjscan.json",
				"${projectDir}/assets/phage_filter_terms_emapper_v2.3.txt",
				params.simple_output
			)

			mge_ch = mge_ch
				.join(mgexpose.out.gff, by: [0, 1], remainder: true)
				.map { speci, genome_id, gdata, flags, mge_gff -> [ speci, genome_id, gdata, flags, mge_gff, null ] }

		} else {
			annotation_data_ch = genomes_ch
				.filter { it[3].PANGENOME_CLUSTERING }
				.map { speci, genome_id, gdata, flags -> [ speci, genome_id, gdata.gff, gdata.conjugation_system_data, gdata.emapper, gdata.gene_clusters, gdata.recombinases, gdata.genome ] }
			annotation_data_ch.dump(pretty: true, tag: "annotation_data_ch")

			mgexpose(
				annotation_data_ch,
				"${projectDir}/assets/mge_rules_ms.txt",
				"${projectDir}/assets/conjscan.json",
				"${projectDir}/assets/phage_filter_terms_emapper_v2.3.txt",
				params.simple_output
			)

			mge_ch = mge_ch
				.join(mgexpose.out.gff, by: [0, 1], remainder: true)
				.join(mgexpose.out.pangenome_info, by: [0, 1], remainder: true)

		} 

		// mge_ch = genomes_ch
		// 	.join(mgexpose.out.gff, by: [0, 1], remainder: true)
		// 	.join(mgexpose.out.pangenome_info, by: [0, 1], remainder: true)
		mge_ch = mge_ch
			.map { speci, genome_id, gdata_old, flags_old, mge_gff, pangenome_info -> 
				def gdata = gdata_old.clone()
				gdata.mge_gff = mge_gff
				gdata.pangenome_info = pangenome_info
				def flags = flags_old.clone()
				flags.MGE_ANNOTATION = (mge_gff != null)
				return [ speci, genome_id, gdata, flags ]
			}

		publish_gene_annotations(
			mge_ch
				.filter { it[2].mge_gff != null }
				.map { speci, genome_id, gdata, flags -> [ speci, genome_id, [ gdata.proteins, gdata.genes, gdata.gff ] ] },
			params.simple_output
		)
	emit:
		genomes = mge_ch
}