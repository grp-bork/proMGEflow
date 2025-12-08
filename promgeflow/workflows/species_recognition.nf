include { recognise_genome; recognise } from "../modules/recognise"

// phasing out nested parameters
params.recognise = [:]
params.recognise.db = null
params.recognise_marker_db = params.recognise.db



workflow species_recognition {

	// run reCOGnise to assign input genome to a specI cluster
	// reCOGnise then automatically collects the specI gene cluster sequences
	// reCOGnise also runs prodigal internally

	take:
		genomes_ch

	main:
		genomes_ch
			.branch {
				annotated: it[1].genes != null
				unannotated: true
			}
			.set{ recognise_input_ch }

		recognise(
			recognise_input_ch.annotated.map { genome_id, gdata -> [ genome_id, gdata.genes, gdata.proteins ] },
			params.recognise_marker_db
		)

		recognise_genome(
			recognise_input_ch.unannotated.map { genome_id, gdata -> [ genome_id, gdata.genome ] },
			params.recognise_marker_db
		)

		genome_speci_ch = recognise_genome.out.genome_speci
			.mix(recognise.out.genome_speci)
			.map { genome_id, file -> [ genome_id, file.text.strip() ] }
			.join(
				recognise_genome.out.speci_status_ok
					.mix(recognise.out.speci_status_ok),
				by: 0,
				remainder: true
			)
			.map { 
				genome_id, speci, status_ok_file -> 
				if (status_ok_file == null) {
					speci = "unknown"
				}
				return [ genome_id, speci ]
			}



		annotations_ch = recognise_genome.out.proteins
			.mix(recognise_genome.out.genes)
			.mix(recognise_genome.out.gff)
			.groupTuple(by: 0, sort: true, size: 3)
			.map { genome_id, annots -> [genome_id, annots[0], annots[1], annots[2] ] }

		recognise_output_ch = recognise_input_ch.unannotated
			.join(annotations_ch, by: 0)
			.map { genome_id, old_gdata, proteins, genes, gff ->
				def gdata = old_gdata.clone()
				gdata.proteins = proteins
				gdata.genes = genes
				gdata.gff = gff
				return [ genome_id, gdata ]
			}
			.mix(recognise_input_ch.annotated)
			.join(genome_speci_ch, by: 0)
			.map { genome_id, gdata_old, speci -> 
				def gdata = gdata_old.clone()
				gdata.speci = speci			
				return [ speci, genome_id, gdata ]
			}
		


			// // .join(genome_speci_ch, by: 0)
			// .join(recognise_input_ch.unannotated, by: 0, )
			// .map { genome_id, annotations, speci -> [speci, genome_id, annotations] }

		// pgenomes_ch = genome_speci_ch
		// 	.join(genomes_ch, by: 0)
		// 	.map { genome_id, speci, genome_fasta -> [speci, genome_id, genome_fasta] }

	emit:
		genomes = recognise_output_ch
		// annotations = annotations_ch
		// genomes = pgenomes_ch

}