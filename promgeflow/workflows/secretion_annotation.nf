include { macsyfinder } from "../modules/macsyfinder"

// phasing out nested parameters
params.txsscan = [:]
params.txsscan.db = null
params.conjscan_models = params.txsscan.db


workflow secretion_annotation {

	take:
		genomes_ch

	main:
		filtered_proteins_ch = genomes_ch.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.proteins ] }

		macsyfinder(filtered_proteins_ch, params.conjscan_models)
		macsy_reports_ch = genomes_ch
			.join(macsyfinder.out.macsy_report, by: [0, 1])
			.map { speci, genome_id, gdata_old, secretion ->
				def gdata = gdata_old.clone()
				gdata.secretion_data = secretion
				return [ speci, genome_id, gdata ]
			}


	emit:
		genomes = macsy_reports_ch

}