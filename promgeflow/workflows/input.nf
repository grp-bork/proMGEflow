params.file_pattern = "**.fna"
def suffix_pattern = params.file_pattern.replaceAll(/\*/, "")

// If the specI classification is known beforehand, it can be passed via the `known_speci` parameter.
// Otherwise, specI will be treated as unknown until it is detected by reCOGnise downstream.
def speci_tag = params.known_speci ?: "unknown"


workflow handle_input_plasmids {

	main:
		genomes_ch = Channel.empty()

		if (params.input_sheet) {
			genomes_ch = Channel
				.fromPath(params.input_sheet)
        		.splitCsv(sep: '\t', header: false)
        		.map { it -> [ "plasmid", it[0], it[1] ] }
						
		} else {
			// Input genomes are genomic fasta files (default ".fna", but can be specified via `file_pattern` parameters) in a directory or directory tree
			// genomes_ch emits tuples (specI, genome_id, genome_fasta)	
			genomes_ch = Channel
				.fromPath(params.input_dir + "/" + params.file_pattern)
				.map { fasta ->
					def genome_id = fasta.name.replaceAll(suffix_pattern, "")
					return tuple("plasmid", genome_id, fasta)
				}
		}

	emit:
		genomes = genomes_ch

}


workflow handle_input_genomes {		

	main:
		genomes_ch = Channel.empty()
		speci_ch = Channel.empty()

		if (params.input_sheet) {
			genomes_ch = Channel
				.fromPath(params.input_sheet)
        		.splitCsv(sep: '\t', header: false)
        		.map { it -> [it[0], it[1], it[2]] }
        
			speci_ch = genomes_ch
				.map {speci, genome, file -> speci}
				.unique()
				.view()
						
		} else {
			// Input genomes are genomic fasta files (default ".fna", but can be specified via `file_pattern` parameters) in a directory or directory tree
			// genomes_ch emits tuples (specI, genome_id, genome_fasta)	
			genomes_ch = Channel
				.fromPath(params.input_dir + "/" + params.file_pattern)
				.map { fasta ->
					def genome_id = fasta.name.replaceAll(suffix_pattern, "")
					return tuple(speci_tag, genome_id, fasta)
				}

			speci_ch = Channel.of(speci_tag)
		}

		genomes_ch
			.branch {
				speci_known: it[0] != "unknown"
				speci_unknown: true
			}
			.set { genomes_speci_ch }


	emit:
		genomes_with_speci = genomes_speci_ch.speci_known
		genomes_without_speci = genomes_speci_ch
			.speci_unknown
			.map { speci, genome_id, genome_fasta -> [genome_id, genome_fasta] }
		speci = speci_ch


}