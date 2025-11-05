include { prodigal; buffered_prodigal } from "../modules/prodigal"

// phasing out prodigal_buffer_size
params.prodigal_buffer_size = -1
params.prodigal_batch_size = params.prodigal_buffer_size


workflow genome_annotation {

	take:
		genomes_ch
	
	main:

		genome_data_ch = genomes_ch.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.genome ] }

	    if (params.prodigal_batch_size != null && params.prodigal_batch_size > 1) {

			def batch_id = 0
			// prodigal_input_ch = genomes_ch
			prodigal_input_ch = genome_data_ch
				.map { speci, genome_id, genome_fasta -> genome_fasta }					
				.buffer(size: params.prodigal_batch_size, remainder: true)
				.map { files -> [ batch_id++, files ] }

			// genome_map = genomes_ch
			genome_map = genome_data_ch 
				.map { speci, genome_id, genome_fasta ->
					// [ genome_fasta.replaceAll(/.+\//, ""), genome_id ]
					[ genome_fasta.getName(), genome_id ]
				}

			genome_map.dump(pretty: true, tag: "genome_map")
			
			prodigal_input_ch.dump(pretty: true, tag: "prodigal_input_ch")

			buffered_prodigal(prodigal_input_ch)			

			buffered_prodigal.out.annotations.dump(pretty: true, tag: "bp_annotations_ch")

			annotations_ch = buffered_prodigal.out.annotations
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
			 		// genomes_ch.map { speci, genome_id, genome_fasta -> [genome_id, speci] },
					genome_data_ch.map { speci, genome_id, genome_fasta -> [genome_id, speci] },
			 		by: 0
				)
			 	.map { genome_id, files, speci -> [speci, genome_id, files] }			

		} else {

			// prodigal(genomes_ch)
			prodigal(genome_data_ch)
			annotations_ch = prodigal.out.proteins
				.mix(prodigal.out.genes)
				.mix(prodigal.out.genome_annotation)
				.groupTuple(by: [0, 1], sort: true, size: 3)
		}

		prodigal_output_ch = genomes_ch
			.join(annotations_ch, by: [0, 1])
			.map { speci, genome_id, gdata_old, proteins, genes, gff ->
				def gdata = gdata_old.clone()
				gdata.genes = genes
				gdata.proteins = proteins
				gdata.gff = gff
				return [ speci, genome_id, gdata ]
			}

	emit:
		// annotations = annotations_ch
		genomes = prodigal_output_ch


}	