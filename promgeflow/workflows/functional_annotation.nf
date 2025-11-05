include { eggnog_mapper } from "../modules/eggnog_mapper"

// phasing out nested parameters and phage-scan-specific parameters
params.phage_scan = [:]
params.phage_scan.db = null
params.emapper_db = params.phage_scan.db


workflow functional_annotation {

	take:
		// filtered_proteins_ch
		genomes_ch

	main:

		filtered_proteins_ch = genomes_ch.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.proteins ] }

		eggnog_mapper(filtered_proteins_ch, params.emapper_db)
		emapper_annotations_ch = eggnog_mapper.out.eggnog
			.join(eggnog_mapper.out.done_sentinel, by: [0, 1])
			.map { speci, genome_id, annotation, sentinel -> [ speci, genome_id, annotation ] }
			.join(genomes_ch, by: [0, 1])
			.map { speci, genome_id, annotation, gdata_old ->
				def gdata = gdata_old.clone()
				gdata.emapper = annotation
				return [ speci, genome_id, gdata ]
			}		

	emit:
		genomes = emapper_annotations_ch	

}