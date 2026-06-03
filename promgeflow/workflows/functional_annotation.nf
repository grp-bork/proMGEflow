include { eggnog_mapper } from "../modules/eggnog_mapper"

// phasing out nested parameters and phage-scan-specific parameters
params.phage_scan = [:]
params.phage_scan.db = null
params.emapper_db = params.phage_scan.db


workflow functional_annotation {

	take:
		genomes_ch

	main:

		filtered_proteins_ch = genomes_ch
			.filter { it[0] != "unknown" && it[2].emapper == null && it[3].RECOMBINASE_SCAN }
			.map { speci, genome_id, gdata, flags -> [ speci, genome_id, gdata.proteins ] }

		eggnog_mapper(filtered_proteins_ch, params.emapper_db)
		emapper_annotations_ch = eggnog_mapper.out.eggnog
			.join(eggnog_mapper.out.done_sentinel, by: [0, 1], remainder: true)
			.map { speci, genome_id, annotation, sentinel -> [ speci, genome_id, (sentinel != null) ? annotation : null ] }
		emapper_annotations_ch = genomes_ch
			.join(emapper_annotations_ch, by: [0, 1], remainder: true)
			.map { speci, genome_id, gdata_old, flags_old, annotation ->
				def gdata = gdata_old.clone()
				if (gdata.emapper == null) {
					gdata.emapper = annotation
				}
				def flags = flags_old.clone()
				flags.FUNCTIONAL_ANNOTATION = (gdata.emapper != null)
				return [ speci, genome_id, gdata, flags ]
			}

	emit:
		genomes = emapper_annotations_ch	

}