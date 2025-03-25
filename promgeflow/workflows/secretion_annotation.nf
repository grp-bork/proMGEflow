include { txsscan } from "../modules/txsscan"


workflow secretion_annotation {

	take:
		filtered_proteins_ch
		genome2speci_map_ch

	main:
		txsscan(filtered_proteins_ch, params.txsscan.db)
		txsscan_reports_ch = txsscan.out.txsscan_report

	emit:
		txsscan = txsscan_reports_ch

}