#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { mgexpose } from "../modules/mgexpose"
include { get_db_seqs } from "../modules/get_db_seqs"

include { genome_annotation } from "./genome_annotation"
include { species_recognition } from "./species_recognition"
include { recombinase_annotation } from "./recombinase_annotation"
include { pangenome_analysis } from "./pangenome_analysis"
include { secretion_annotation } from "./secretion_annotation"
include { functional_annotation } from "./functional_annotation"

include { handle_input_genomes } from "./input"


params.file_pattern = "**.fna"
params.genome_buffer_size = 100
print "PARAMS:\n" + params


// def suffix_pattern = params.file_pattern.replaceAll(/\*/, "")

// // If the specI classification is known beforehand, it can be passed via the `known_speci` parameter.
// // Otherwise, specI will be treated as unknown until it is detected by reCOGnise downstream.
// def speci_tag = params.known_speci ?: "unknown_speci"

workflow full_annotation {

	classified_genomes_ch = Channel.empty()
	speci_seqs_ch = Channel.empty()
	speci_source_ch = Channel.empty()

	// prodigal output channels
	// pproteins_ch = Channel.empty()
	// pgenes_ch = Channel.empty()
	// pgffs_ch = Channel.empty()

	handle_input_genomes()

	// // Input genomes are genomic fasta files (default ".fna", but can be specified via `file_pattern` parameters) in a directory or directory tree
	// // genomes_ch emits tuples (specI, genome_id, genome_fasta)	
	// genomes_ch = Channel
	// 	.fromPath(params.input_dir + "/" + params.file_pattern)
	// 	.map { fasta ->
	// 		def genome_id = fasta.name.replaceAll(suffix_pattern, "")
	// 		return tuple(speci_tag, genome_id, fasta)
	// 	}

	/* STEP 1A: genome annotation via prodigal if speci is known */
	genome_annotation(handle_input_genomes.out.genomes_with_speci)
	/* STEP 1B: specI assignment via reCOGnise */
	species_recognition(handle_input_genomes.out.genomes_without_speci)

	genomes_ch = handle_input_genomes.out.genomes_with_speci
		.mix(species_recognition.out.genomes)
	genomes_ch.dump(pretty: true, tag: "genomes_ch")

	// we don't need this, we derive it from filtered genomes below
	// speci_ch = handle_input_genomes.out.speci
	// 	.mix(species_recognition.out.speci)
	// 	.unique()

	// prodigal output channels
	// pproteins_ch = genome_annotation.out.proteins
	// 	.mix(species_recognition.out.proteins)
	// pgenes_ch = genome_annotation.out.genes
	// 	.mix(species_recognition.out.genes)
	// pgffs_ch = genome_annotation.out.gffs
	// 	.mix(species_recognition.out.gffs)	

	annotations_ch = genome_annotation.out.annotations
		.mix(species_recognition.out.annotations)

	annotations_ch.dump(pretty: true, tag: "annotations_ch")
	
	// if (params.known_speci) {

	// 	/* STEP 1A: genome annotation via prodigal if speci is known */
	// 	genome_annotation(genomes_ch)

	// 	pproteins_ch = genome_annotation.out.proteins
	// 	pgenes_ch = genome_annotation.out.genes
	// 	pgffs_ch = genome_annotation.out.gffs

	// } else {

	// 	/* STEP 1B: specI assignment via reCOGnise */

	// 	// run reCOGnise to assign input genome to a specI cluster
	// 	// reCOGnise then automatically collects the specI gene cluster sequences
	// 	// reCOGnise also runs prodigal internally
	// 	species_recognition(genomes_ch.map { speci, genome_id, genome_fasta -> [genome_id, genome_fasta] })

	// 	pproteins_ch = species_recognition.out.proteins
	// 	pgenes_ch = species_recognition.out.genes
	// 	pgffs_ch = species_recognition.out.gffs
	// 	genomes_ch = species_recognition.out.genomes

	// }
	// pproteins_ch.dump(pretty: true, tag: "pproteins_ch")
	
	// genome2speci_map_ch = pgenes_ch
	// 	.map { speci, genome_id, genes -> return tuple(genome_id, speci) }

	// genome2speci_map_ch.dump(pretty: true, tag: "genome2speci_map_ch")

	/* STEP 2: Run recombinase annotation */
	recombinase_annotation(
		annotations_ch
			.map { it -> [it[0], it[1], it[2][0]] }
	)
	// recombinase_annotation(pproteins_ch) //, genome2speci_map_ch)
	recombinase_annotation.out.recombinases.dump(pretty: true, tag: "recombinases")

	filtered_ch = annotations_ch
		.join(recombinase_annotation.out.recombinases, by: [0, 1])
		.map { speci, genome_id, annotations, recombinases -> [speci, genome_id, annotations] }
	
	filtered_ch.dump(pretty: true, tag: "filtered_ch")

	speci_seqs_ch = filtered_ch
		.map { speci, genome_id, annotations -> speci }
		.filter { it != "unknown" }
		.unique()
	speci_seqs_ch.dump(pretty: true, tag: "speci_seqs_ch")
	
	// filtered_genes_ch = pgenes_ch
	// 	.join(recombinase_annotation.out.recombinases, by: [0, 1])
	// 		.map { speci, genome_id, genes, recombinases ->
	// 			return tuple(speci, genome_id, genes) 
	// 		}
	// filtered_proteins_ch = pproteins_ch
	// 	.join(recombinase_annotation.out.recombinases, by: [0, 1])
	// 		.map { speci, genome_id, proteins, recombinases ->
	// 			return tuple(speci, genome_id, proteins) 
	// 		}
	// filtered_gff_ch = pgffs_ch
	// 	.join(recombinase_annotation.out.recombinases, by: [0, 1])
	// 		.map { speci, genome_id, gff, recombinases -> [speci, genome_id, gff] }

	// filtered_genes_ch.dump(pretty: true, tag: "filtered_genes_ch")

	// speci_seqs_ch = filtered_genes_ch
	//  	.map { speci, genome_id, genes -> return speci }
	// 	.filter { it != "unknown" }  // should not happen, but let's be defensive.
	// 	.unique()

	// speci_seqs_ch.dump(pretty: true, tag: "speci_seqs_ch")
	
	// get_db_seqs(speci_seqs_ch, "progenomes3_db", params.genedata.db, params.genedata.db_credentials, params.genedata.cache)
	// params.gene_cluster_seqdb = "/g/bork6/schudoma/experiments/mge_refseqindex/sp095_refdb/sp095_refdb.tar"
	get_db_seqs(speci_seqs_ch, params.gene_cluster_seqdb)

	filtered_proteins_ch = filtered_ch
		.map { speci, genome_id, annotations -> [speci, genome_id, annotations[0]] }
	filtered_genes_ch = filtered_ch
		.map { speci, genome_id, annotations -> [speci, genome_id, annotations[1]] }
	filtered_gff_ch = filtered_ch
		.map { speci, genome_id, annotations -> [speci, genome_id, annotations[2]] }

	/* STEP 3 Perform gene clustering */
	pangenome_analysis(filtered_genes_ch, get_db_seqs.out.sequences)

	/* STEP 4 Protein annotation - phage signals and secretion systems */
	secretion_annotation(filtered_proteins_ch) //, genome2speci_map_ch)
	functional_annotation(filtered_proteins_ch) //, genome2speci_map_ch)
	
	/* STEP 5 Annotate the genomes with island data and assign mges */
	
	// [speci, bin_id, gene_coords, txsscan, emapper, clusters, recombinases, genome_fa]
	annotation_data_ch = filtered_gff_ch
		.join( secretion_annotation.out.txsscan, by: [0, 1] )
		.join( functional_annotation.out.annotation, by: [0, 1] )
		.join( pangenome_analysis.out.clusters, by: [0, 1] )
		.join( recombinase_annotation.out.recombinases, by: [0, 1] )
		.join( genomes_ch, by: [0, 1] )

	mgexpose(
		annotation_data_ch,
		"${projectDir}/assets/mge_rules_ms.txt",
		"${projectDir}/assets/txsscan_rules.txt",
		"${projectDir}/assets/phage_filter_terms_emapper_v2.3.txt"
	)

}
