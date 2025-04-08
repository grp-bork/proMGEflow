include { recognise_genome } from "../modules/recognise"


workflow species_recognition {

	// run reCOGnise to assign input genome to a specI cluster
	// reCOGnise then automatically collects the specI gene cluster sequences
	// reCOGnise also runs prodigal internally

	take:
		genomes_ch

	main:
		recognise_genome(
			genomes_ch,
			params.recognise.db
		)

		genome_speci_ch = recognise_genome.out.speci_status_ok
			.join(
				recognise_genome.out.genome_speci
					.map { genome_id, file -> return tuple(genome_id, file.text.strip()) },
				by: 0
			)
			.map { genome_id, file, speci -> return tuple(genome_id, speci) }

		annotations_ch = recognise_genome.out.proteins
			.mix(recognise_genome.out.genes)
			.mix(recognise_genome.out.gff)
			.groupTuple(by: [0, 1], sort: true)
			.join(genome_speci_ch, by: 0)
			.map { genome_id, annotations, speci -> [speci, genome_id, annotations] }

		pgenomes_ch = genome_speci_ch
			.join(genomes_ch, by: 0)
			.map { genome_id, speci, genome_fasta -> [speci, genome_id, genome_fasta] }

	emit:
		annotations = annotations_ch
		genomes = pgenomes_ch

}