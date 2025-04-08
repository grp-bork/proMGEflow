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
				.map { annotation_file -> 
					[ annotation_file.getName().replaceAll(/\.(faa|ffn|gff)$/, ""), annotation_file ]
				}
				.groupTuple(by: 0, sort: true)
				.join(
					genomes_ch.map { speci, genome_id, genome_fasta -> [genome_id, speci] },
					by: 0
				)
				.map { genome_id, annotation_file, speci -> [speci, genome_id, annotation_file] }			

		} else {

			prodigal(genomes_ch)
			annotations_ch = prodigal.out.proteins
				.mix(prodigal.out.genes)
				.mix(prodigal.out.genome_annotation)
				.groupTuple(by: [0, 1], sort: true)
		}

	emit:
		annotations = annotations_ch


}