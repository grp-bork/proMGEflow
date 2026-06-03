include { macsyfinder } from "../modules/macsyfinder"

// phasing out nested parameters
params.txsscan = [:]
params.txsscan.db = params.txsscan_db
params.conjscan_models = params.txsscan.db


workflow conjugation_system_annotation {

	take:
		genomes_ch

	main:
		filtered_proteins_ch = genomes_ch
			// in certain situations, we want to annotate conjugation systems
			// on genomes without pangenome
			// e.g. bulk annotation with preliminary pangenome data
			.filter { it[3].PANGENOME_ESTIMATION || params.force_conjugation_system_analysis }
			.map { speci, genome_id, gdata, flags -> [ speci, genome_id, gdata.proteins ] }

		macsyfinder(filtered_proteins_ch, params.conjscan_models)
		macsy_reports_ch = genomes_ch
			.join(macsyfinder.out.macsy_report, by: [0, 1], remainder: true)
			.map { speci, genome_id, gdata_old, flags_old, conjugation_system ->
				def gdata = gdata_old.clone()
				gdata.conjugation_system_data = conjugation_system
				if (!flags_old.PANGENOME_ESTIMATION) {
					// little hack for no-pangenome runs
					gdata.gene_clusters = file("${workDir}/DUMMY_CLUSTERS.txt")
				}
				def flags = flags_old.clone()
				flags.CONJUGATION_SYSTEM_ANNOTATION = (gdata.conjugation_system_data != null && gdata.conjugation_system_data.text.strip() != "")
				return [ speci, genome_id, gdata, flags ]
			}

	emit:
		genomes = macsy_reports_ch

}