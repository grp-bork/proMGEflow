include { recombinase_scan } from "../modules/recombinase_scan"
// include { publish_recombinase_scan } from "../modules/publish"

// phasing out nested parameters
params.recombinase_scan = [:]
params.recombinase_scan.db = null
params.recombinase_scan_db = params.recombinase_scan.db


workflow recombinase_annotation {

	take:
		genomes_ch

	main:
		genomes_ch.dump(pretty: true, tag: "recombinase_annotation_input")

		proteins_ch = genomes_ch
			.map { speci, genome_id, gdata, flags -> [ speci, genome_id, gdata.proteins ] }
			.filter { it[2] != null }

		proteins_ch.dump(pretty: true, tag: "proteins_ch_in_recombinase_annotation")

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
		recombinase_output_ch = genomes_ch
			.join(recombinase_output_ch, by: [0, 1], remainder: true)
			.map { speci, genome_id, gdata_old, flags_old, recombinases, recomb_table, recomb_gff ->
				def gdata = gdata_old.clone()
				gdata.recombinases = recombinases
				gdata.recomb_table = recomb_table
				gdata.recomb_gff = recomb_gff
				def flags = flags_old.clone()
				flags.RECOMBINASE_SCAN = (recombinases != null && recomb_table != null && recomb_gff != null)
				return [ speci, genome_id, gdata, flags ]
			}

		// publish_recombinase_scan(
		// 	recombinase_output_ch
		// 		.filter { it[2].recomb_table != null && it[2].recomb_gff != null }
		// 		.map { speci, genome_id, gdata, flags -> [ speci, genome_id, gdata.recomb_table, gdata.recomb_gff ] },
		// 	params.simple_output
		// )

	emit:
		genomes = recombinase_output_ch

}
