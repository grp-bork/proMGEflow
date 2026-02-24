#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { mgexpose } from "../modules/mgexpose"
include { get_db_seqs } from "../modules/get_db_seqs"
include { publish_gene_annotations; publish_recombinase_scan } from "../modules/publish"

include { genome_annotation } from "./genome_annotation"
include { species_recognition } from "./species_recognition"
include { recombinase_annotation } from "./recombinase_annotation"
include { pangenome_analysis } from "./pangenome_analysis"
include { secretion_annotation; secretion_annotation as forced_secretion_annotation } from "./secretion_annotation"
include { functional_annotation } from "./functional_annotation"

include { handle_input_genomes } from "./input"


params.genome_buffer_size = 100
print "PARAMS:\n" + params


workflow guided_annotation {

	handle_input_genomes()

	

	/* genome annotation via prodigal */
	// genome_annotation(
	// 	handle_input_genomes.out.to_species_recognition.mix(handle_input_genomes.out.to_genome_annotation)
	// )
	
	// // prodigal output channels
	// genomes_ch = genome_annotation.out.genomes
	// 	// .mix(handle_input_genomes.out.to_recombinase_scan)
	// genomes_ch.dump(pretty: true, tag: "genomes_ch")

	// /* STEP 2: Run recombinase annotation */
	// recombinase_annotation(genomes_ch)

	// /* STEP 2b: Filter by recombinase presence */
	// with_recombinase_ch = recombinase_annotation.out.genomes



	/* STEP Y Publish recombinase annotations */

	// publish_recombinase_scan(
	// 	recombinase_annotation.out.genomes
	// 		.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.recomb_table, gdata.recomb_gff ] },
	// 	params.simple_output
	// )

}
