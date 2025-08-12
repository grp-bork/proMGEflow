include { recombinase_scan } from "../modules/recombinase_scan"

// phasing out nested parameters
params.recombinase_scan = [:]
params.recombinase_scan.db = null
params.recombinase_scan_db = params.recombinase_scan.db


workflow recombinase_annotation {

	take:
		pproteins_ch

	main:		
		recombinase_scan(
			pproteins_ch,
			params.recombinase_scan_db,
			"${projectDir}/assets/mge_rules_ms.txt"
		)
		annotated_recombinases_ch = recombinase_scan.out.recombinases
		mge_predictions_ch = recombinase_scan.out.recomb_table
		mge_predictions_gff_ch = recombinase_scan.out.mge_pred_gff

	emit:
		recombinases = annotated_recombinases_ch
		mge_predictions = mge_predictions_ch
		mge_predictions_gff = mge_predictions_gff_ch
}
