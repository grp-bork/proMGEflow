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


process extract_recombinase_contigs {
	container "quay.io/biocontainers/seqtk:1.5--h577a1d6_1"
	memory { 4.GB * task.attempt }
	time { 2.h * task.attempt }
	tag "${genome_id}"

	input:
	tuple val(genome_id), path(fasta), path(gff)

	output:
	tuple val(genome_id), path("${genome_id}.recombinase_contigs.fa.gz"), emit: contigs

	script:
	"""
	seqtk subseq ${fasta} <(grep -v "^#" ${gff} | cut -f 1 | uniq | sort -u) | gzip -c - > ${genome_id}.recombinase_contigs.fa.gz
	"""
	// seqtk subseq ${sample.id}_1.fastq chimeras.txt >> chimeras.fastq

}

params.mgedb = "/g/scb2/bork/grekova/projects/promge_website/14112025/dataset/mge100/mge_sequences_unique.fa"
params.minimap_x = "asm20"

process map_mgedb {
	container "quay.io/biocontainers/minimap2:2.30--h577a1d6_0"
	// container "registry.git.embl.org/schudoma/align-docker:latest"
	memory { 32.GB * task.attempt }
	time { 8.h * task.attempt }
	cpus 4
	tag "${genome_id}"


	input:
	tuple val(genome_id), path(fasta)
	tuple val(db_id), path(db)

	script:
	"""
	minimap2 -x ${params.minimap_x} -t ${task.cpus} -a -c -L --eqx --sam-hit-only -o ${genome_id}.${db_id}.sam ${db} ${fasta}
	"""


}


workflow guided_annotation {

	handle_input_genomes()

	// handle_input_genomes().out.to_genome_annotation.dump(pretty: true, tag: "tga_ch")


	/* genome annotation via prodigal */
	genome_annotation(handle_input_genomes.out.to_species_recognition.map { genome_id, gdata -> [ "unknown", genome_id, gdata ] } )
	// genome_annotation(
	// 	handle_input_genomes.out.to_species_recognition.mix(handle_input_genomes.out.to_genome_annotation)
	// )
	
	// // prodigal output channels
	genomes_ch = genome_annotation.out.genomes
	// 	// .mix(handle_input_genomes.out.to_recombinase_scan)
	genomes_ch.dump(pretty: true, tag: "genomes_ch")

	/* STEP 2: Run recombinase annotation */
	recombinase_annotation(genomes_ch)

	/* STEP 2b: Filter by recombinase presence */
	with_recombinase_ch = recombinase_annotation.out.genomes

	extract_recombinase_contigs(
		with_recombinase_ch
			.map { speci, genome_id, gdata -> [ genome_id, gdata.genome, gdata.recomb_gff ] }
	)

	map_mgedb(
		extract_recombinase_contigs.out.contigs,
		["promge100", params.mgedb]
	)


	/* STEP Y Publish recombinase annotations */

	// publish_recombinase_scan(
	// 	recombinase_annotation.out.genomes
	// 		.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.recomb_table, gdata.recomb_gff ] },
	// 	params.simple_output
	// )

}
