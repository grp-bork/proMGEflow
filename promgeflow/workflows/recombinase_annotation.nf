include { recombinase_scan } from "../modules/recombinase_scan"


workflow recombinase_annotation {

	take:
		pproteins_ch
		// genome2speci_map_ch

	main:		
		recombinase_scan(
			pproteins_ch,
			params.recombinase_scan.db,
			"${projectDir}/assets/mge_rules_ms.txt"
		)
		annotated_recombinases_ch = recombinase_scan.out.recombinases

	emit:
		recombinases = annotated_recombinases_ch

}