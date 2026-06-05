#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { mgexpose_denovo } from "./mgexpose"
include { get_db_seqs } from "../modules/get_db_seqs"

include { genome_annotation } from "./genome_annotation"
include { species_recognition } from "./species_recognition"
include { recombinase_annotation } from "./recombinase_annotation"
include { pangenome_analysis } from "./pangenome_analysis"
include { conjugation_system_annotation; conjugation_system_annotation as forced_conjugation_system_annotation } from "./conjugation_system_annotation"
include { functional_annotation } from "./functional_annotation"
include { summarise_and_publish } from "./summarise"

include { handle_input_genomes } from "./input"


params.genome_buffer_size = 100
print "PARAMS:\n" + params


workflow full_annotation {

	handle_input_genomes()

	/* STEP 1A: genome annotation via prodigal for genomes with known speci */
	genome_annotation(handle_input_genomes.out.to_genome_annotation)
	
	/* STEP 1B: specI assignment via reCOGnise otherwise */
	species_recognition(handle_input_genomes.out.to_species_recognition)

	// prodigal output channels
	genomes_ch = genome_annotation.out.genomes
		.mix(species_recognition.out.genomes)
		.mix(handle_input_genomes.out.to_recombinase_scan)
	genomes_ch.dump(pretty: true, tag: "genomes_ch")

	/* STEP 2: Run recombinase annotation */
	recombinase_annotation(genomes_ch)

	functional_annotation(recombinase_annotation.out.genomes)

	/* STEP 2c: Obtain speci reference gene sequences */
	speci_seqs_ch = functional_annotation.out.genomes
		.map { speci, genome_id, gdata, flags -> speci }
		.filter { it != "unknown" }
		.unique()
	
	// params.gene_cluster_seqdb = "/g/bork6/schudoma/experiments/mge_refseqindex/sp095_refdb/sp095_refdb.tar"
	get_db_seqs(speci_seqs_ch, params.gene_cluster_seqdb)
	speci_refseqs_ch = get_db_seqs.out.sequences
		.join(get_db_seqs.out.done_sentinel, by: 0)
		.map { speci, sequences, sentinel -> [ speci, sequences ] }

	with_functional_annotation_ch = functional_annotation.out.genomes
		.join(
			functional_annotation.out.genomes
				.map { speci, genome_id, gdata, flags -> [ speci, genome_id ] }
				.combine(speci_refseqs_ch, by: 0),
			by: [0, 1], remainder: true
		)
		.map { speci, genome_id, gdata, flags_old, sequences -> 
			def flags = flags_old.clone()
			flags.SPECI_CLUSTER_SEQS = (sequences != null)
			return [ speci, genome_id, gdata, flags ]
		}

	/* STEP 3 Perform gene clustering */
	pangenome_analysis(
		with_functional_annotation_ch,
		speci_refseqs_ch
	)

	/* STEP 4 Protein annotation - phage signals and conjugation systems */
	conjugation_system_annotation(pangenome_analysis.out.genomes)
	
	/* STEP 5 Annotate the genomes with island data and assign mges */
	mgexpose_denovo(conjugation_system_annotation.out.genomes)

	summarise_and_publish(mgexpose_denovo.out.genomes)

	// if (true) {

	// 	pub_recombinases_ch = recombinase_annotation.out.mge_predictions
	// 		.mix(recombinase_annotation.out.mge_predictions_gff)
	// 		.groupTuple(by: [0, 1], size: 2)

	// 	pub_recombinases_nospeci_ch = pub_recombinases_ch
	// 		.filter { it[0] == "unknown" }

	// 	pub_recombinases_nomge_ch = pub_recombinases_ch
	// 		.filter { it[0] != "unknown" }
	// 		.join(mgexpose.out.gff, by: [0, 1], remainder: true)
	// 		.filter { it[3] == null }			
	// 		.map { speci, genome_id, recombinases, no_mge -> [speci, genome_id, recombinases ] }

	// 	publish_ch = annotations_ch
	// 		.join(
	// 			mgexpose.out.gff.mix(mgexpose.out.fasta).groupTuple(by: [0, 1], size: 2),
	// 			by: [0, 1]
	// 		)
	// 		.map { speci, genome_id, annotations, mges -> [ speci, genome_id, [ annotations[0], annotations[1], annotations[2], mges[0], mges[1] ] ] }
	// 		.mix(pub_recombinases_nospeci_ch)
	// 		.mix(pub_recombinases_nomge_ch)

	// 	publish_ch.dump(pretty: true, tag: "publish_ch")

	// 	publish_results(publish_ch, params.simple_output, params.tarball_output)

	// }


}
