include { prodigal; buffered_prodigal; pyrodigal; buffered_pyrodigal } from "../modules/prodigal"

params.prodigal_buffer_size = -1


workflow genome_annotation {

	take:
		genomes_ch
	
	main:

	    if (params.prodigal_buffer_size != null && params.prodigal_buffer_size > 1) {

			def batch_id = 0
			prodigal_input_ch = genomes_ch
				.map { speci, genome_id, genome_fasta -> genome_fasta }					
				.buffer(size: params.prodigal_buffer_size, remainder: true)
				.map { files -> [ batch_id++, files ] }

			genome_map = genomes_ch
				.map { speci, genome_id, genome_fasta ->
					[ genome_fasta.replaceAll(/.+\//, ""), genome_id ]
				}

			genome_map.dump(pretty: true, tag: "genome_map")
			
			prodigal_input_ch.dump(pretty: true, tag: "prodigal_input_ch")

			buffered_pyrodigal(prodigal_input_ch)			

			buffered_pyrodigal.out.annotations.dump(pretty: true, tag: "bp_annotations_ch")

			annotations_ch = buffered_pyrodigal.out.annotations
				.flatten()
				.map { annotation_file -> 
					[ annotation_file.getName().replaceAll(/\.(faa|ffn|gff)$/, ""), annotation_file ]
				}
				.groupTuple(by: 0, sort: true, size: 3)
				.join(genome_map, by: 0)
				.map { fn, files, genome_id -> [ genome_id, files ] }

			annotations_ch.dump(pretty: true, tag: "annotations_post_ch")
			annotations_ch = annotations_ch
				.join(
			 		genomes_ch.map { speci, genome_id, genome_fasta -> [genome_id, speci] },
			 		by: 0
				)
			 	.map { genome_id, files, speci -> [speci, genome_id, files] }			

		} else {

			pyrodigal(genomes_ch)
			annotations_ch = pyrodigal.out.proteins
				.mix(pyrodigal.out.genes)
				.mix(pyrodigal.out.genome_annotation)
				.groupTuple(by: [0, 1], sort: true, size: 3)
		}

	emit:
		annotations = annotations_ch


}	