include { txsscan } from "../modules/txsscan"

// phasing out nested parameters
params.txsscan ?: [:]
params.txsscan.db = null
params.txsscan_db = params.txsscan.db


workflow secretion_annotation {

	take:
		filtered_proteins_ch

	main:
		txsscan(filtered_proteins_ch, params.txsscan_db)
		txsscan_reports_ch = txsscan.out.txsscan_report

	emit:
		txsscan = txsscan_reports_ch

}