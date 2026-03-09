#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
	tuple val(genome_id), path("${genome_id}.recombinase_contigs.fa.gz"), emit: contigs, optional: true
	tuple val(genome_id), path("${genome_id}.RECOMBINASE_CONTIGS.DONE"), emit: sentinel

	script:
	"""
	seqtk subseq ${fasta} <(grep -v "^#" ${gff} | cut -f 1 | uniq | sort -u) | gzip -c - > ${genome_id}.recombinase_contigs.fa.gz

	if [[ -z \$(zcat ${genome_id}.recombinase_contigs.fa.gz | head -n 1) ]]; then rm -fv ${genome_id}.recombinase_contigs.fa.gz; fi

	touch ${genome_id}.RECOMBINASE_CONTIGS.DONE
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
	tuple val(genome_id), path("${genome_id}.${db_id}.sam.gz"), emit: sam
	tuple val(genome_id), path("${genome_id}.mmi"), emit: index

	script:
	"""
	minimap2 -x ${params.minimap_x} -d ${genome_id}.mmi ${fasta}

	minimap2 -x ${params.minimap_x} -t ${task.cpus} -a -c -L --eqx --sam-hit-only ${genome_id}.mmi ${db} | samtools sort -O SAM -o ${genome_id}.${db_id}.samx
	awk -v OFS='\t' '/^[^@]/ {\$10="*"; print \$0}' ${genome_id}.${db_id}.samx | gzip -c - > ${genome_id}.${db_id}.sam.gz
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
	guide_extract_matches.py ${sam} > \$(basename ${sam} .sam.gz).bed
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

process mgexpose {
	label "annotate_genome"
	label "medium"
	// container "ghcr.io/cschu/mgexpose:v3.8.0"
	tag "${genome_id}"

	input:
	tuple val(genome_id), path(gff), path(recombinases), path(genome_fa), path(emapper), path(conjscan)
	path(mge_rules)
	path(conjscan_rules)
	path(phage_filter_terms)

	output:
	tuple val(genome_id), path("**/*.mge_islands.gff3"), emit: gff, optional: true
	path("**/*.NO_MGE"), emit: no_mge, optional: true
	path("**/*.mge_islands.ffn.gz"), emit: fasta, optional: true
	path("**/*.gene_info.txt"), emit: gene_info, optional: true
	
	script:
	def outdir = "${genome_id}"

	"""
	mkdir -p ${outdir}/

	if [[ "${gff}" == *".gz" ]]; then
		gzip -dc ${gff} > mgexpose.gff
	else
		ln -sf ${gff} mgexpose.gff
	fi


	echo mgexpose reannotate ${genome_id} mgexpose.gff \
			--annotation_mode raw_islands \
			--recombinase_hits ${recombinases} \
			--mge_rules ${mge_rules} \
			--txs_macsy_rules ${conjscan_rules} \
			--txs_macsy_report ${conjscan} \
			--phage_eggnog_data ${emapper} \
			--phage_filter_terms ${phage_filter_terms} \
			--output_dir ${outdir} \
			--extract_islands ${genome_fa} \
			--output_suffix mge_islands
	mgexpose reannotate ${genome_id} mgexpose.gff \
			--annotation_mode raw_islands \
			--recombinase_hits ${recombinases} \
			--mge_rules ${mge_rules} \
			--txs_macsy_rules ${conjscan_rules} \
			--txs_macsy_report ${conjscan} \
			--phage_eggnog_data ${emapper} \
			--phage_filter_terms ${phage_filter_terms} \
			--output_dir ${outdir} \
			--extract_islands ${genome_fa} \
			--output_suffix mge_islands 

	islands_gff=${outdir}/${genome_id}.mge_islands.gff3
	(grep mobile_genetic_element \${islands_gff} | grep -v mge= > ${genome_id}.NO_MGE) || true
	if [[ -s ${genome_id}.NO_MGE ]]; then 
		mv -v ${genome_id}.NO_MGE ${outdir}/
	else
		rm -f ${genome_id}.NO_MGE
	fi

	rm -vf mgexpose.gff
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

	secretion_annotation(downstream_ch)

	mgexpose_ch = convert_to_gff_and_extract_proteins.out.gff
		.join(
			recombinase_annotation.out.genomes
				.map { speci, genome_id, gdata -> [ genome_id, gdata.recombinases, gdata.genome ] },
			by: 0
		)
		.join(
			functional_annotation.out.genomes
				.map { speci, genome_id, gdata -> [ genome_id, gdata.emapper ] },
			by: 0
		)
		.join(
			secretion_annotation.out.genomes
				.map { speci, genome_id, gdata -> [ genome_id, gdata.secretion_data, ] },
				by: 0
		)

	// tuple val(genome_id), path(gff), path(recombinases), path(genome_fa), path(emapper), path(conjscan)

	mgexpose(
		mgexpose_ch,
		"${projectDir}/assets/mge_rules_ms.txt",
		"${projectDir}/assets/conjscan.json",
		"${projectDir}/assets/phage_filter_terms_emapper_v2.3.txt"
	)

	
	// mgexpose reannotate -o . --annotation_mode raw_islands 
	// --recombinase_hits work/83/067fefaefa515a9538f6bbf25f7d37/unknown/SAMN10093300-assembled/SAMN10093300-assembled.recombinase_hmmsearch.besthits.out
	// --mge_rules ~/.nextflow/assets/grp-bork/promgeflow/assets/mge_rules_ms.txt 
	// --txs_macsy_report work/a8/48a2b4f5dbe4817696dbc06997ffa5/unknown/SAMN10093300-assembled/SAMN10093300-assembled.all_systems.tsv 
	// --txs_macsy_rules ~/.nextflow/assets/grp-bork/promgeflow/assets/conjscan.json 
	// --phage_eggnog_data work/94/f2664c75b2236a507282b4e84ddabe/unknown/SAMN10093300-assembled/SAMN10093300-assembled.emapper.annotations 
	// --phage_filter_terms ~/.nextflow/assets/grp-bork/promgeflow/assets/phage_filter_terms_emapper_v2.3.txt 
	// --extract_islands input/SAMN10093300-assembled.fa.gz 
	// work/b7/48f9c5ef89b96909cbc51ba93fad72/SAMN10093300-assembled.mge_candidates.gff3 GENOME2 > log.txt
	


	/* STEP Y Publish recombinase annotations */

	// publish_recombinase_scan(
	// 	recombinase_annotation.out.genomes
	// 		.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.recomb_table, gdata.recomb_gff ] },
	// 	params.simple_output
	// )

}
