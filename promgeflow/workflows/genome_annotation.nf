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

			def batch_id = 0
			buffered_prodigal(
				genomes_ch
					.map { speci, genome_id, genome_fasta -> genome_fasta }
					.buffer(size: params.prodigal_buffer_size, remainder: true)
					.map { files -> [ batch_id++, files ] },
				suffix_pattern
			)
			annotations_ch = buffered_prodigal.out.annotations
				.flatten()
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