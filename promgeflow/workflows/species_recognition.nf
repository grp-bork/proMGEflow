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

		pproteins_ch = recognise_genome.out.proteins
			.join( genome_speci_ch, by: 0)
			.map { genome_id, proteins, speci -> return tuple(speci, genome_id, proteins) }
		pproteins_ch.dump(pretty: true, tag: "pproteins_ch")

		pgenes_ch = recognise_genome.out.genes
			.join( genome_speci_ch, by: 0)
			.map { genome_id, genes, speci -> return tuple(speci, genome_id, genes) }
		
		pgffs_ch = recognise_genome.out.gff
			.join( genome_speci_ch, by: 0)
			.map { genome_id, gff, speci -> return tuple(speci, genome_id, gff) }

		pgenomes_ch = genomes_ch
			.join( genomes_ch, by: 0 )
			.map { genome_id, speci, genome_fasta -> [speci, genome_id, genome_fasta] }
	
		mixed_ch = pproteins_ch
				.mix(pgenes_ch)
				.mix(pgffs_ch)
				.groupTuple(by: [0, 1], sort: true)


	emit:
		proteins = pproteins_ch
		genes = pgenes_ch
		gffs = pgffs_ch
		genomes = pgenomes_ch
		speci = pgenomes_ch
			.map { speci, genome_id, genome_fasta -> speci }
			.unique()
		mixed = mixed_ch

}