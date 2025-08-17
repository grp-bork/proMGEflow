#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { mgexpose } from "../modules/mgexpose"
include { get_db_seqs } from "../modules/get_db_seqs"
include { publish_gene_annotations; publish_recombinase_scan; publish_results } from "../modules/publish"

include { genome_annotation } from "./genome_annotation"
include { species_recognition } from "./species_recognition"
include { recombinase_annotation } from "./recombinase_annotation"
include { pangenome_analysis } from "./pangenome_analysis"
include { secretion_annotation } from "./secretion_annotation"
include { functional_annotation } from "./functional_annotation"

include { handle_input_genomes } from "./input"


params.genome_buffer_size = 100
print "PARAMS:\n" + params


workflow full_annotation {

	handle_input_genomes()

	/* STEP 1A: genome annotation via prodigal for genomes with known speci */
	genome_annotation(handle_input_genomes.out.genomes_with_speci)
	/* STEP 1B: specI assignment via reCOGnise otherwise */
	species_recognition(handle_input_genomes.out.genomes_without_speci)

	genomes_ch = handle_input_genomes.out.genomes_with_speci
		.mix(species_recognition.out.genomes)
	genomes_ch.dump(pretty: true, tag: "genomes_ch")

	// prodigal output channels
	annotations_ch = genome_annotation.out.annotations
		.mix(species_recognition.out.annotations)

	annotations_ch.dump(pretty: true, tag: "annotations_ch")

	/* STEP 2: Run recombinase annotation */
	recombinase_annotation(
		annotations_ch
			.map { it -> [it[0], it[1], it[2][0]] }
	)

	/* STEP 2b: Filter by recombinase presence */
	with_recombinase_ch = annotations_ch
		.join(recombinase_annotation.out.recombinases, by: [0, 1])
		.filter { it[0] != "unknown" }
		.map { speci, genome_id, annotations, recombinases -> [speci, genome_id, annotations] }

	functional_annotation(
		with_recombinase_ch.map { speci, genome_id, annotations -> [speci, genome_id, annotations[0]] }
	)

	with_functional_annotation_ch = with_recombinase_ch
		.join(functional_annotation.out.annotation, by: [0, 1])
		.map { speci, genome_id, annotations, functional_annotation -> [ speci, genome_id, annotations ] }

	/* STEP 2c: Obtain speci reference gene sequences */
	speci_seqs_ch = with_functional_annotation_ch
		.map { speci, genome_id, annotations -> speci }
		.filter { it != "unknown" }
		.unique()
	
	// params.gene_cluster_seqdb = "/g/bork6/schudoma/experiments/mge_refseqindex/sp095_refdb/sp095_refdb.tar"
	get_db_seqs(speci_seqs_ch, params.gene_cluster_seqdb)
	speci_refseqs_ch = get_db_seqs.out.sequences
		.join(get_db_seqs.out.done_sentinel, by: 0)
		.map { speci, sequences, sentinel -> [ speci, sequences ] }

	/* STEP 3 Perform gene clustering */
	pangenome_analysis(
		with_functional_annotation_ch.map { speci, genome_id, annotations -> [ speci, genome_id, annotations[1] ] },
		speci_refseqs_ch
	)

	with_cluster_ch = with_functional_annotation_ch
		.join(pangenome_analysis.out.clusters, by: [0, 1])
		.map { speci, genome_id, annotations, clusters -> [ speci, genome_id, annotations ] }

	/* STEP 4 Protein annotation - phage signals and secretion systems */
	secretion_annotation(with_cluster_ch.map { speci, genome_id, annotations -> [ speci, genome_id, annotations[0] ] })
	
	/* STEP 5 Annotate the genomes with island data and assign mges */
	annotation_data_ch = with_cluster_ch.map { speci, genome_id, annotations -> [ speci, genome_id, annotations[2] ] }
		.join( secretion_annotation.out.txsscan, by: [0, 1] )
		.join( functional_annotation.out.annotation, by: [0, 1] )
		.join( pangenome_analysis.out.clusters, by: [0, 1] )
		.join( recombinase_annotation.out.recombinases, by: [0, 1] )
		.join( genomes_ch, by: [0, 1] )

	mgexpose(
		annotation_data_ch,
		"${projectDir}/assets/mge_rules_ms.txt",
		"${projectDir}/assets/txsscan_rules.txt",
		"${projectDir}/assets/phage_filter_terms_emapper_v2.3.txt",
		params.simple_output
	)

	if (true) {

		pub_recombinases_ch = recombinase_annotation.out.mge_predictions
			.mix(recombinase_annotation.out.mge_predictions_gff)
			.filter { it[0] == "unknown" }
			.groupTuple(by: [0, 1], size: 2)
		// pub_recombinases_ch.dump(pretty: true, tag: "pub_recombinases_ch")
		// [DUMP: pub_recombinases_ch] [
 		//    "unknown",
    	// 	"NT12426",
		// 	[
		// 		"/g/bork6/schudoma/projects/mge/typas_hogan_TEST_DONT_TOUCH/work/f7/c18fae0ed1453f5856c80f41560c59/unknown/NT12426/NT12426.recombinase_based_MGE_predictions.tsv",
		// 		"/g/bork6/schudoma/projects/mge/typas_hogan_TEST_DONT_TOUCH/work/f7/c18fae0ed1453f5856c80f41560c59/unknown/NT12426/NT12426.predicted_recombinase_mges.gff3"
		// 	]
		// ]

		publish_ch = annotations_ch
			.join(
				mgexpose.out.gff.mix(mgexpose.out.fasta).groupTuple(by: [0, 1], size: 2),
				by: [0, 1]
			)
			.map { speci, genome_id, annotations, mges -> [ speci, genome_id, [ annotations[0], annotations[1], annotations[2], mges[0], mges[1] ] ] }
			.mix(pub_recombinases_ch)

		// [DUMP: publish_ch] [
		// 	"specI_v4_00001",
		// 	"NT12014",
		// 	[
		// 		"/g/bork6/schudoma/projects/mge/typas_hogan_TEST_DONT_TOUCH/work/50/ab817e0603271f117769bd8e24d6be/specI_v4_00001/NT12014/NT12014.faa",
		// 		"/g/bork6/schudoma/projects/mge/typas_hogan_TEST_DONT_TOUCH/work/50/ab817e0603271f117769bd8e24d6be/specI_v4_00001/NT12014/NT12014.ffn",
		// 		"/g/bork6/schudoma/projects/mge/typas_hogan_TEST_DONT_TOUCH/work/50/ab817e0603271f117769bd8e24d6be/specI_v4_00001/NT12014/NT12014.gff"
		// 	],
		// 	[
		// 		"/g/bork6/schudoma/projects/mge/typas_hogan_TEST_DONT_TOUCH/work/87/58a40b5d8f2fcdd755bc8ca9173045/NT12014/NT12014.mge_islands.gff3",
		// 		"/g/bork6/schudoma/projects/mge/typas_hogan_TEST_DONT_TOUCH/work/87/58a40b5d8f2fcdd755bc8ca9173045/NT12014/NT12014.mge_islands.ffn.gz"
    	// ]

					// 
			// .mix
			// .groupTuple(
			// 	recombinase_annotation.out.mge_predictions
			// 		.join(recombinase_annotation.out.mge_predictions_gff, by: [0, 1])
			// 		.filter { it[0] == "unknown" }
			// )
		publish_ch.dump(pretty: true, tag: "publish_ch")

		// publish_results(publish_ch, params.simple_output, params.tarball_output)

	} else {

		publish_gene_annotations(
			annotations_ch
				.join(mgexpose.out.gff, by: [0, 1])
				.map { speci, genome_id, annotations, mge_gff -> [ speci, genome_id, annotations ] },
			params.simple_output
		)

		publish_recombinase_scan(
			recombinase_annotation.out.mge_predictions
				.join(recombinase_annotation.out.mge_predictions_gff, by: [0, 1])
				.filter { it[0] == "unknown" },
			params.simple_output
		)

	}


}
