#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { mgexpose_region } from "../modules/mgexpose"

include { mgexpose_denovo } from "./mgexpose"
include { genome_annotation } from "./genome_annotation"
include { recombinase_annotation } from "./recombinase_annotation"
include { conjugation_system_annotation } from "./conjugation_system_annotation"
include { functional_annotation } from "./functional_annotation"
include { summarise_and_publish } from "./summarise"
include { handle_input_contigs } from "./input"


params.genome_buffer_size = 100
print "PARAMS:\n" + params


workflow contig_annotation {

	handle_input_contigs()

	/* STEP 1A: genome annotation via prodigal for genomes with known speci */
	genome_annotation(handle_input_contigs.out.genomes)	

	// prodigal output channels
	genomes_ch = genome_annotation.out.genomes

	genomes_ch.dump(pretty: true, tag: "genomes_ch")

	/* STEP 2: Run recombinase annotation */
	recombinase_annotation(genomes_ch)

	functional_annotation(recombinase_annotation.out.genomes)

	conjugation_system_annotation(functional_annotation.out.genomes)
	
	/* STEP 5 Annotate the genomes with island data and assign mges */
	mgexpose_denovo(conjugation_system_annotation.out.genomes)

	summarise_and_publish(mgexpose_denovo.out.genomes)

	// annotation_data_ch = conjugation_system_annotation.out.genomes
	// 	.map { speci, genome_id, gdata -> [ speci, genome_id, "dummy_region", gdata.gff, gdata.conjugation_system_data, gdata.emapper, gdata.recombinases, gdata.genome ] }

	// annotation_data_ch.dump(pretty: true, tag: "annotation_data_ch")

	// mgexpose_region(
	// 	annotation_data_ch,
	// 	"${projectDir}/assets/mge_rules_ms.txt",
	// 	"${projectDir}/assets/conjscan.json",
	// 	"${projectDir}/assets/phage_filter_terms_emapper_v2.3.txt",
	// 	params.simple_output
	// )

}
