#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { mgexpose_region } from "../modules/mgexpose"
include { publish_gene_annotations; publish_recombinase_scan } from "../modules/publish"

include { genome_annotation } from "./genome_annotation"
include { recombinase_annotation } from "./recombinase_annotation"
include { secretion_annotation; secretion_annotation as forced_secretion_annotation } from "./secretion_annotation"
include { functional_annotation } from "./functional_annotation"

include { handle_input_plasmids } from "./input"


params.genome_buffer_size = 100
print "PARAMS:\n" + params


workflow plasmid_annotation {

	handle_input_plasmids()

	genomes_ch = handle_input_plasmids.out.genomes

	genomes_ch.dump(pretty: true, tag: "genomes_ch")

	/* STEP 1A: genome annotation via prodigal for genomes with known speci */
	genome_annotation(genomes_ch)	

	// prodigal output channels
		annotations_ch = genome_annotation.out.genomes

	annotations_ch.dump(pretty: true, tag: "annotations_ch")

	/* STEP 2: Run recombinase annotation */
	recombinase_annotation(
		annotations_ch
		//	.map { it -> [it[0], it[1], it[2][0]] }
	)

	with_recombinase_ch = recombinase_annotation.out.genomes

	with_recombinase_ch
		.branch {
			has_emapper: it[2].emapper != null
			to_emapper: true
		}
		.set { emapper_input_ch }

	functional_annotation(emapper_input_ch.to_emapper)

	with_functional_annotation_ch = emapper_input_ch.has_emapper
		.mix(functional_annotation.out.genomes)

	secretion_annotation(with_functional_annotation_ch)
	secretion_ch = secretion_annotation.out.genomes
	if (params.force_secretion_analysis) {
		// in certain situations, we want to annotate the secretion system 
		// on genomes without pangenome
		// e.g. bulk annotation with preliminary pangenome data
		dummy_clusters = file('DUMMY_CLUSTERS.txt')


		forced_secretion_annotation(with_functional_annotation_ch)
		forced_secretion_ch = forced_secretion_annotation.out.genomes
			.map { speci, genome_id, gdata_old -> 
				def gdata = gdata_old.clone()
				gdata.gene_clusters = dummy_clusters
				return [ speci, genome_id, gdata ]
			}
		secretion_ch = secretion_ch.mix(forced_secretion_ch)
	}

	/* STEP 5 Annotate the genomes with island data and assign mges */
	// 	tuple val(speci), val(genome_id), path(gff), path(txsscan), path(emapper), path(gene_clusters), path(recombinases), path(genome_fa)

	annotation_data_ch = secretion_ch
		.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.region_id, gdata.gff, gdata.secretion_data, gdata.emapper, gdata.recombinases, gdata.genome ] }

	// annotation_data_ch = filtered_gff_ch
	// 	.join( secretion_annotation.out.txsscan, by: [0, 1] )
	// 	.join( functional_annotation.out.annotation, by: [0, 1] )
	// 	.join( handle_input_plasmids.out.regions, by: [0, 1] )
	// 	.join( recombinase_annotation.out.recombinases, by: [0, 1] )
	// 	.join( genomes_ch, by: [0, 1] )

	annotation_data_ch.dump(pretty: true, tag: "annotation_data_ch")

	mgexpose_region(
		annotation_data_ch,
		"${projectDir}/assets/mge_rules_ms.txt",
		"${projectDir}/assets/txsscan_rules.txt",
		"${projectDir}/assets/phage_filter_terms_emapper_v2.3.txt",
		params.simple_output
	)

	/* STEP X Publish gene annotations of input genomes that were not pre-annotated */

	publish_gene_annotations(
		secretion_annotation.out.genomes
			.join(mgexpose_region.out.gff, by: [0, 1])
			.join(genomes_ch, by: [0, 1])
			.map { speci, genome_id, gdata, mge_gff, gdata_raw -> [ speci, genome_id, [ gdata.proteins, gdata.genes, gdata.gff ] ] },
		params.simple_output
	)

	/* STEP Y Publish recombinase annotations */

	publish_recombinase_scan(
		recombinase_annotation.out.genomes
			.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.recomb_table, gdata.recomb_gff ] },
		params.simple_output
	)

}
