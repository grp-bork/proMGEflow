#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { mgexpose_region } from "../modules/mgexpose"
include { get_db_seqs } from "../modules/get_db_seqs"
include { publish_gene_annotations } from "../modules/publish"

include { genome_annotation } from "./genome_annotation"
include { recombinase_annotation } from "./recombinase_annotation"
include { pangenome_analysis } from "./pangenome_analysis"
include { secretion_annotation } from "./secretion_annotation"
include { functional_annotation } from "./functional_annotation"

include { handle_input_plasmids } from "./input"


params.genome_buffer_size = 100
print "PARAMS:\n" + params


workflow plasmid_annotation {

	handle_input_plasmids()

	genomes_ch = handle_input_plasmids.out.genomes

	genomes_ch.dump(pretty: true, tag: "genomes_ch")

	/* STEP 1A: genome annotation via prodigal for genomes with known speci */
	genome_annotation(handle_input_plasmids.out.genomes)	

	// prodigal output channels
	annotations_ch = genome_annotation.out.annotations

	// publish_annotations(annotations_ch)

	annotations_ch.dump(pretty: true, tag: "annotations_ch")

	/* STEP 2: Run recombinase annotation */
	recombinase_annotation(
		annotations_ch
			.map { it -> [it[0], it[1], it[2][0]] }
	)

	/* STEP 2b: Filter by recombinase presence */
	filtered_ch = annotations_ch
		.join(recombinase_annotation.out.recombinases, by: [0, 1])
		.map { speci, genome_id, annotations, recombinases -> [speci, genome_id, annotations] }

	filtered_proteins_ch = filtered_ch
		.map { speci, genome_id, annotations -> [speci, genome_id, annotations[0]] }
	filtered_genes_ch = filtered_ch
		.map { speci, genome_id, annotations -> [speci, genome_id, annotations[1]] }
	filtered_gff_ch = filtered_ch
		.map { speci, genome_id, annotations -> [speci, genome_id, annotations[2]] }

	secretion_annotation(filtered_proteins_ch)
	functional_annotation(filtered_proteins_ch)
	
	/* STEP 5 Annotate the genomes with island data and assign mges */
	// 	tuple val(speci), val(genome_id), path(gff), path(txsscan), path(emapper), path(gene_clusters), path(recombinases), path(genome_fa)

	annotation_data_ch = filtered_gff_ch
		.join( secretion_annotation.out.txsscan, by: [0, 1] )
		.join( functional_annotation.out.annotation, by: [0, 1] )
		.join( handle_input_plasmids.out.regions, by: [0, 1] )
		.join( recombinase_annotation.out.recombinases, by: [0, 1] )
		.join( genomes_ch, by: [0, 1] )

	annotation_data_ch.dump(pretty: true, tag: "annotation_data_ch")

	mgexpose_region(
		annotation_data_ch,
		"${projectDir}/assets/mge_rules_ms.txt",
		"${projectDir}/assets/txsscan_rules.txt",
		"${projectDir}/assets/phage_filter_terms_emapper_v2.3.txt"
	)

	mgexpose_region.out.gff.dump(pretty: true, tag: "mgexpose_region_ch")

	publish_ch = annotations_ch
			.join(mgexpose_region.out.gff, by: [0, 1])
			.map { speci, genome_id, annotations, mge_gff -> [ speci, genome_id, annotations ] }

	publish_ch.dump(pretty: true, tag: "final_annotations_ch")

	publish_annotations(
		publish_ch
	)

}
