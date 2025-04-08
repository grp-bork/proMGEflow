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
			annotations_ch = buffered_prodigal.out.annotations
				.flatten()
				.map { genome_fasta -> 
					[ genome_fasta.getName().replaceAll(/\.(faa|ffn|gff)$/, ""), genome_fasta ]
				}
				.join(
					genomes_ch.map { speci, genome_id, genome_fasta -> [genome_id, speci] },
					by: 0
				)
				.map { genome_id, genome_fasta, speci -> [speci, genome_id, genome_fasta] }

				// .map { file -> [
				// 	params.known_speci, file.getName().replaceAll(/\.(faa|ffn|gff)$/, ""), file
				// ]}
			pproteins_ch = annotations_ch
				.filter { it[2].getName().endsWith(".faa") }
			pgenes_ch = annotations_ch
				.filter { it[2].getName().endsWith(".ffn") }
			pgffs_ch = annotations_ch
				.filter { it[2].getName().endsWith(".gff") }

		} else {

			prodigal(genomes_ch)
			pproteins_ch = prodigal.out.proteins
			pgenes_ch = prodigal.out.genes
			pgffs_ch = prodigal.out.genome_annotation

		}

		mixed_ch = pproteins_ch
				.mix(pgenes_ch)
				.mix(pgffs_ch)
				.groupTuple(by: [0, 1], sort: true)

	emit:
		proteins = pproteins_ch
		genes = pgenes_ch
		gffs = pgffs_ch
		mixed = mixed_ch


}