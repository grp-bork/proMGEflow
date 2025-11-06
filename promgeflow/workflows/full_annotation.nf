#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { mgexpose } from "../modules/mgexpose"
include { get_db_seqs } from "../modules/get_db_seqs"
include { publish_gene_annotations; publish_recombinase_scan } from "../modules/publish"

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
	// genome_annotation(handle_input_genomes.out.genomes_with_speci)
	genome_annotation(handle_input_genomes.out.to_genome_annotation)
	/* STEP 1B: specI assignment via reCOGnise otherwise */
	// species_recognition(handle_input_genomes.out.genomes_without_speci)
	species_recognition(handle_input_genomes.out.to_species_recognition)

	// genomes_ch = handle_input_genomes.out.genomes_with_speci
	// 	.mix(species_recognition.out.genomes)
	// potentially no longer needed
	xgenomes_ch = handle_input_genomes.out.to_genome_annotation
		.mix(handle_input_genomes.out.to_species_recognition)
		.mix(handle_input_genomes.out.to_recombinase_scan)
	xgenomes_ch.dump(pretty: true, tag: "xgenomes_ch")

	// prodigal output channels
	genomes_ch = genome_annotation.out.genomes
		.mix(species_recognition.out.genomes)
		.mix(handle_input_genomes.out.to_recombinase_scan)
	genomes_ch.dump(pretty: true, tag: "genomes_ch")
	// annotations_ch.dump(pretty: true, tag: "annotations_ch")

	/* STEP 2: Run recombinase annotation */
	recombinase_annotation(genomes_ch)
		
		// annotations_ch
		// 	.map { it -> [it[0], it[1], it[2][0]] }
	

	/* STEP 2b: Filter by recombinase presence */
	// with_recombinase_ch = annotations_ch
	with_speci_and_recombinase_ch = recombinase_annotation.out.genomes
		// .join(recombinase_annotation.out.recombinases, by: [0, 1])
		.filter { it[0] != "unknown" }
		// .map { speci, genome_id, annotations, recombinases -> [speci, genome_id, annotations] }

	with_speci_and_recombinase_ch
		.branch {
			has_emapper: it[2].emapper != null
			to_emapper: true
		}
		.set { emapper_input_ch }

	functional_annotation(emapper_input_ch.to_emapper)
		// with_recombinase_ch.map { speci, genome_id, annotations -> [speci, genome_id, annotations[0]] }
		// with_speci_and_recombinase_ch

	with_functional_annotation_ch = emapper_input_ch.has_emapper
		.mix(functional_annotation.out.genomes)

	// with_functional_annotation_ch = with_recombinase_ch
	// 	.join(functional_annotation.out.annotation, by: [0, 1])
	// 	.map { speci, genome_id, annotations, functional_annotation -> [ speci, genome_id, annotations ] }

	/* STEP 2c: Obtain speci reference gene sequences */
	speci_seqs_ch = with_functional_annotation_ch
		// .map { speci, genome_id, annotations -> speci }
		.map { speci, genome_id, gdata -> speci }
		.filter { it != "unknown" }
		.unique()
	
	// params.gene_cluster_seqdb = "/g/bork6/schudoma/experiments/mge_refseqindex/sp095_refdb/sp095_refdb.tar"
	get_db_seqs(speci_seqs_ch, params.gene_cluster_seqdb)
	speci_refseqs_ch = get_db_seqs.out.sequences
		.join(get_db_seqs.out.done_sentinel, by: 0)
		.map { speci, sequences, sentinel -> [ speci, sequences ] }

	/* STEP 3 Perform gene clustering */
	pangenome_analysis(
		// with_functional_annotation_ch.map { speci, genome_id, annotations -> [ speci, genome_id, annotations[1] ] },
		with_functional_annotation_ch,
		speci_refseqs_ch
	)

	// with_cluster_ch = with_functional_annotation_ch
	// 	.join(pangenome_analysis.out.clusters, by: [0, 1])
	// 	.map { speci, genome_id, annotations, clusters -> [ speci, genome_id, annotations ] }
	// with_clusters_ch = pangenome_analysis.out.genomes

	/* STEP 4 Protein annotation - phage signals and secretion systems */
	// secretion_annotation(with_cluster_ch.map { speci, genome_id, annotations -> [ speci, genome_id, annotations[0] ] })
	secretion_annotation(pangenome_analysis.out.genomes)
	
	/* STEP 5 Annotate the genomes with island data and assign mges */
	// tuple val(speci), val(genome_id), path(gff), path(txsscan), path(emapper), path(gene_clusters), path(recombinases), path(genome_fa)

	annotation_data_ch = secretion_annotation.out.genomes
		.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.gff, gdata.secretion, gdata.emapper, gdata.gene_clusters, gdata.recombinases, gdata.genome ] }

	annotation_data_ch.dump(pretty: true, tag: "annotation_data_ch")
	// annotation_data_ch = with_cluster_ch.map { speci, genome_id, annotations -> [ speci, genome_id, annotations[2] ] }
	// 	.join( secretion_annotation.out.txsscan, by: [0, 1] )
	// 	.join( functional_annotation.out.annotation, by: [0, 1] )
	// 	.join( pangenome_analysis.out.clusters, by: [0, 1] )
	// 	.join( recombinase_annotation.out.recombinases, by: [0, 1] )
	// 	.join( genomes_ch, by: [0, 1] )

	mgexpose(
		annotation_data_ch,
		"${projectDir}/assets/mge_rules_ms.txt",
		"${projectDir}/assets/txsscan_rules.txt",
		"${projectDir}/assets/phage_filter_terms_emapper_v2.3.txt",
		params.simple_output
	)

	publish_gene_annotations(
		// annotations_ch
		secretion_annotation.out.genomes
			.join(mgexpose.out.gff, by: [0, 1])
			// .map { speci, genome_id, annotations, mge_gff -> [ speci, genome_id, annotations ] },
			.map { speci, genome_id, gdata -> [ speci, genome_id, [ gdata.proteins, gdata.genes, gdata.gff ] ] },
		params.simple_output
	)

	publish_recombinase_scan(
		with_speci_and_recombinase_ch.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.recomb_table, gdata.recomb_gff ] },
		// recombinase_annotation.out.mge_predictions
		// 	.join(recombinase_annotation.out.mge_predictions_gff, by: [0, 1])
		// 	.filter { it[0] == "unknown" },
		params.simple_output
	)

}
