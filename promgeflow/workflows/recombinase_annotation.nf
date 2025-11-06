include { recombinase_scan } from "../modules/recombinase_scan"

// phasing out nested parameters
params.recombinase_scan = [:]
params.recombinase_scan.db = null
params.recombinase_scan_db = params.recombinase_scan.db


workflow recombinase_annotation {

	take:
		// pproteins_ch
		genomes_ch

	main:

		proteins_ch = genomes_ch
			.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.proteins ] }

		recombinase_scan(
			proteins_ch,
			params.recombinase_scan_db,
			"${projectDir}/assets/mge_rules_ms.txt"
		)

		recombinase_output_ch = recombinase_scan.out.done_sentinel
			.join(recombinase_scan.out.recombinases, by: [0, 1])
			.join(recombinase_scan.out.recomb_table, by: [0, 1])
			.join(recombinase_scan.out.mge_pred_gff, by: [0, 1])
			.map { speci, genome_id, sentinel, recombinases, recomb_table, recomb_gff -> 
				return [ speci, genome_id, recombinases, recomb_table, recomb_gff ]
			}
			.join(genomes_ch, by: [0, 1])
			.map { speci, genome_id, sentinel, recombinases, recomb_table, recomb_gff, old_gdata ->
				def gdata = old_gdata.clone()
				gdata.recombinases = recombinases
				gdata.recomb_table = recomb_table
				gdata.recomb_gff = recomb_gff
				return [ speci, genome_id, gdata ]
			}
		
	emit:
		genomes = recombinase_output_ch
		// recombinases = annotated_recombinases_ch
		// mge_predictions = mge_predictions_ch
		// mge_predictions_gff = mge_predictions_gff_ch
}
