include { txsscan } from "../modules/txsscan"

// phasing out nested parameters
params.txsscan = [:]
params.txsscan.db = null
params.txsscan_db = params.txsscan.db


workflow secretion_annotation {

	take:
		// filtered_proteins_ch
		genomes_ch

	main:
		filtered_proteins_ch = genomes_ch.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.proteins ] }

		txsscan(filtered_proteins_ch, params.txsscan_db)
		txsscan_reports_ch = genomes_ch
			.join(txsscan.out.txsscan_report, by: [0, 1])
			.map { speci, genome_id, gdata_old, secretion ->
				def gdata = gdata_old.clone()
				gdata.secretion_data = secretion
				return [ speci, genome_id, gdata ]
			}


	emit:
		genomes = txsscan_reports_ch

}