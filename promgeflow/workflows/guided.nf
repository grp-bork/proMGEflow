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

params.mgedb = "/scratch/schudoma/databases/mge/mge_sequences_unique.fa"
params.minimap_x = "asm20"

process map_mgedb {
	// container "quay.io/biocontainers/minimap2:2.30--h577a1d6_0"
	container "registry.git.embl.org/schudoma/align-docker:with_minimap2"
	memory { 32.GB * task.attempt }
	time { 8.h * task.attempt }
	cpus 4
	tag "${genome_id}"


	input:
	tuple val(genome_id), path(fasta)
	tuple val(db_id), path(db)

	output:
	tuple val(genome_id), path("${genome_id}.${db_id}.sam"), emit: sam
	tuple val(genome_id), path("${genome_id}.mmi"), emit: index

	script:
	"""
	minimap2 -x ${params.minimap_x} -d ${genome_id}.mmi ${fasta}

	minimap2 -x ${params.minimap_x} -t ${task.cpus} -a -c -L --eqx --sam-hit-only ${genome_id}.mmi ${db} | samtools sort -O SAM -o ${genome_id}.${db_id}.sam 
	"""


}

process extract_matches {
	time { 30.m * task.attempt }
	tag "${genome_id}"

	input:
	tuple val(genome_id), path(sam)

	output:
	tuple val(genome_id), path("*.bed"), emit: bed

	script:
	"""
	guide_extract_matches.py ${sam} > \$(basename ${sam} .sam).bed
	"""


}

process check_recombinase_hits {
	container "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
	time { 30.m * task.attempt }
	tag "${genome_id}"

	input:
	tuple val(genome_id), path(bed), path(gff)

	output:
	tuple val(genome_id), path("${genome_id}.recombinase_hits.tsv"), emit: results

	script:
	"""
	bedtools intersect -wo -a ${bed} -b ${gff} > ${genome_id}.recombinase_hits.tsv
	"""
}

process extract_mge_candidates {
	container "ghcr.io/grp-bork/mgexpose:v3.8.0"
	time { 2.h * task.attempt }
	memory { 4.GB * task.attempt }
	tag "${genome_id}"

	input:
	tuple val(genome_id), path(table)

	output:
	tuple val(genome_id), path("${genome_id}.mge_candidates.bed"), emit: bed

	script:
	"""
	sort -k1,1 -k8,8g -k9,9g ${table} | guide.py - ${genome_id}
	"""
}

process add_genes {
	container "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
	time { 2.h * task.attempt }
	memory { 4.GB * task.attempt }
	tag "${genome_id}"

	input:
	tuple val(genome_id), path(bed), path(gff)

	output:
	tuple val(genome_id), path("${genome_id}.mge_candidates.with_genes.tsv"), emit: table

	script:
	"""
	head -1 ${gff} > tmp.gff
	grep -v "^#" ${gff} >> tmp.gff

	bedtools intersect -wo -a ${bed} -b tmp.gff | sort -k1,1h -k2,2h -k3,3h > ${genome_id}.mge_candidates.with_genes.tsv

	rm -frv tmp.gff
	"""
}

process convert_to_gff_and_extract_proteins {
	container "ghcr.io/grp-bork/mgexpose:v3.8.0"
	time { 2.h * task.attempt }
	memory { 8.GB * task.attempt }
	tag "${genome_id}"
	// maxForks 1

	input:
	tuple val(genome_id), path(table), path(faa)

	output:
	tuple val(genome_id), path("${genome_id}.cargo.faa"), emit: proteins
	tuple val(genome_id), path("${genome_id}.mge_candidates.gff3"), emit: gff

	script:
	"""
	guide_extract_regions.py ${table} ${faa} ${genome_id}
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

	extract_matches(map_mgedb.out.sam)


	ch = extract_matches.out.bed.join(
		with_recombinase_ch
			.map { speci, genome_id, gdata -> [ genome_id, gdata.recomb_gff ] },
		by: 0
	) 

	check_recombinase_hits(ch)

	extract_mge_candidates(check_recombinase_hits.out.results)


	add_genes(
		extract_mge_candidates.out.bed
			.join(
				with_recombinase_ch.map { speci, genome_id, gdata -> [ genome_id, gdata.gff ] },
				by: 0
			)
	)

	convert_to_gff_and_extract_proteins(
		add_genes.out.table
			.join(
				with_recombinase_ch.map { speci, genome_id, gdata -> [ genome_id, gdata.proteins ] },
				by: 0
			)
	)

	downstream_ch = convert_to_gff_and_extract_proteins.out.proteins
		.map { genome_id, proteins ->
			def gdata = [:]
			gdata.proteins = proteins
			return [ "unknown", genome_id, gdata ]
		}

		// .join(with_recombinase_ch.map { speci, genome_id, gdata -> [ genome_id, gdata ] }, by: 0)
		// .map {
		// 	genome_id, proteins, old_gdata ->
		// 		def gdata = old_gdata.clone()
		// 		gdata.proteins = proteins
		// 		return [ "unknown", genome_id, gdata ]
		// 	}
			

	functional_annotation(downstream_ch)

	// secretion_annotation(

	// )


	/* STEP Y Publish recombinase annotations */

	// publish_recombinase_scan(
	// 	recombinase_annotation.out.genomes
	// 		.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.recomb_table, gdata.recomb_gff ] },
	// 	params.simple_output
	// )

}
