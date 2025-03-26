include { prodigal; buffered_prodigal } from "../modules/prodigal"

params.prodigal_buffer_size = -1


workflow genome_annotation {

	take:
		genomes_ch
	
	main:
	    // suffix_pattern and speci_tag are defined in main.nf but not recognized otherwise
	    def suffix_pattern = params.file_pattern.replaceAll(/\*/, "")
	    def speci_tag = params.known_speci ?: "unknown_speci"

		if (params.prodigal_buffer_size != null && params.prodigal_buffer_size > 1) {

			buffered_prodigal(genomes_ch.buffer(size: params.prodigal_buffer_size, remainder: true))
			annotations_ch = buffered_prodigal.out.annotations
				.map { file -> [
					file.getName().replaceAll(/\.(faa|ffn|gff)$/, ""), file
				]}
			pproteins_ch = annotations_ch
				.filter { it[1].getName().endsWith(".faa") }
			pgenes_ch = annotations_ch
				.filter { it[1].getName().endsWith(".ffn") }
			pgffs_ch = annotations_ch
				.filter { it[1].getName().endsWith(".gff") }

		} else {

			prodigal(genomes_ch)
			pproteins_ch = prodigal.out.proteins
			pgenes_ch = prodigal.out.genes
			pgffs_ch = prodigal.out.genome_annotation

		}

	emit:
		proteins = pproteins_ch
		genes = pgenes_ch
		gffs = pgffs_ch

}