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


workflow full_annotation {

	handle_input_genomes()

	/* STEP 1A: genome annotation via prodigal for genomes with known speci */
	genome_annotation(handle_input_genomes.out.to_genome_annotation)
	/* STEP 1B: specI assignment via reCOGnise otherwise */
	species_recognition(handle_input_genomes.out.to_species_recognition)

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

	/* STEP 2: Run recombinase annotation */
	recombinase_annotation(genomes_ch)


	/* STEP 2b: Filter by recombinase presence */
	with_speci_and_recombinase_ch = recombinase_annotation.out.genomes
		.filter { it[0] != "unknown" }

	with_speci_and_recombinase_ch
		.branch {
			has_emapper: it[2].emapper != null
			to_emapper: true
		}
		.set { emapper_input_ch }

	functional_annotation(emapper_input_ch.to_emapper)

	with_functional_annotation_ch = emapper_input_ch.has_emapper
		.mix(functional_annotation.out.genomes)

	/* STEP 2c: Obtain speci reference gene sequences */
	speci_seqs_ch = with_functional_annotation_ch
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
		with_functional_annotation_ch,
		speci_refseqs_ch
	)

	/* STEP 4 Protein annotation - phage signals and secretion systems */
	// secretion_annotation(with_cluster_ch.map { speci, genome_id, annotations -> [ speci, genome_id, annotations[0] ] })
	secretion_annotation(pangenome_analysis.out.genomes)
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
	// tuple val(speci), val(genome_id), path(gff), path(secretion), path(emapper), path(gene_clusters), path(recombinases), path(genome_fa)

	annotation_data_ch = secretion_ch
		.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.gff, gdata.secretion_data, gdata.emapper, gdata.gene_clusters, gdata.recombinases, gdata.genome ] }

	annotation_data_ch.dump(pretty: true, tag: "annotation_data_ch")

	mgexpose(
		annotation_data_ch,
		"${projectDir}/assets/mge_rules_ms.txt",
		"${projectDir}/assets/conjscan.json",
		"${projectDir}/assets/phage_filter_terms_emapper_v2.3.txt",
		params.simple_output
	)

	pangenome_ch = Channel.fromPath("${projectDir}/assets/speci_sizes_pg3.txt")
		.splitCsv(header: false, sep: '\t')
		.join(
			mgexpose.out.pangenome_info.splitCsv(header: false, sep: '\t'),
			by: 0
		)
		.map { speci, n_genomes, genome, data -> [ speci, genome, n_genomes, data[2], data[3], data[4], ((data[3].toFloat() / data[2].toFloat()) * 100.0).round(2) ]}
	pangenome_ch.dump(pretty: true, tag: "pangenome_ch")

	[DUMP: pangenome_ch] [
    "specI_v4_08282",
    "2",
    "GCA_000014125.1",
    [
        "specI_v4_08282",
        "GCA_000014125.1",
        "2321",
        "711",
        "1609"
    ]
]
	
	// pangenome_ch.collectFile(name: "pangenome_info.txt", newLine: true, storeDir: params.output_dir)


	publish_gene_annotations(
		secretion_annotation.out.genomes
			.join(mgexpose.out.gff, by: [0, 1])
			.join(handle_input_genomes.out.to_genome_annotation, by: [0, 1])
			.map { speci, genome_id, gdata, mge_gff, gdata_raw -> [ speci, genome_id, [ gdata.proteins, gdata.genes, gdata.gff ] ] },
		params.simple_output
	)

	publish_recombinase_scan(
		recombinase_annotation.out.genomes
			.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.recomb_table, gdata.recomb_gff ] },
		params.simple_output
	)

}
