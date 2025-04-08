include { prodigal; buffered_prodigal } from "../modules/prodigal"

params.prodigal_buffer_size = -1


workflow genome_annotation {

	take:
		genomes_ch
	
	main:
	    def suffix_pattern = params.file_pattern.replaceAll(/\*/, "")
	    
		if (params.prodigal_buffer_size != null && params.prodigal_buffer_size > 1) {

			def batch_id = 0
			buffered_prodigal(
				genomes_ch
					.map { speci, genome_id, genome_fasta -> genome_fasta }
					.buffer(size: params.prodigal_buffer_size, remainder: true)
					.map { files -> [ batch_id++, files ] },
				suffix_pattern
			)
			buffered_prodigal.out.annotations.dump(pretty: true, tag: "buffered_prodigal")
			annotations_ch = buffered_prodigal.out.annotations
				.flatten()
				.map { annotation_file -> 
					[ annotation_file.getName().replaceAll(/\.(faa|ffn|gff)$/, ""), annotation_file ]
				}
				.join(
					genomes_ch.map { speci, genome_id, genome_fasta -> [genome_id, speci] },
					by: 0
				)
				.map { genome_id, annotation_file, speci -> [speci, genome_id, annotation_file] }


			annotations_ch.dump(pretty: true, tag: "annotations_ch")

				// .map { file -> [
				// 	params.known_speci, file.getName().replaceAll(/\.(faa|ffn|gff)$/, ""), file
				// ]}
			pproteins_ch = annotations_ch
				.filter { it[2].getName().endsWith(".faa") }
			pgenes_ch = annotations_ch
				.filter { it[2].getName().endsWith(".ffn") }
			pgffs_ch = annotations_ch
				.filter { it[2].getName().endsWith(".gff") }
			mixed_ch = annotations_ch
				.groupTuple(by: [0, 1], sort: true)
			mixed_ch.dump(pretty: true, tag: "grouped_annotations_ch")

		} else {

			prodigal(genomes_ch)
			pproteins_ch = prodigal.out.proteins
			pgenes_ch = prodigal.out.genes
			pgffs_ch = prodigal.out.genome_annotation
			mixed_ch = Channel.empty() //TODO

		}

		// mixed_ch = pproteins_ch
		// 		.join(pgenes_ch, by: [0, 1])
		// 		.join(pgffs_ch, by: [0, 1])
		// 		// .groupTuple(by: [0, 1], sort: true)

	emit:
		proteins = pproteins_ch
		genes = pgenes_ch
		gffs = pgffs_ch
		mixed = mixed_ch


}