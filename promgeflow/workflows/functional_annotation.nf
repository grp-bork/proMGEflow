include { eggnog_mapper } from "../modules/eggnog_mapper"


workflow functional_annotation {

	take:
		filtered_proteins_ch
		genome2speci_map_ch

	main:
		eggnog_mapper(filtered_proteins_ch, params.phage_scan.db)
		emapper_annotations_ch = eggnog_mapper.out.eggnog

	emit:
		annotation = emapper_annotations_ch	

}