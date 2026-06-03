// If the specI classification is known beforehand, it can be passed via the `known_speci` parameter.
// Otherwise, specI will be treated as unknown until it is detected by reCOGnise downstream.
def speci_tag = params.known_speci ?: "unknown"


workflow handle_input_contigs {

	main:
		genomes_ch = Channel.empty()

		if (params.input_sheet) {
			genomes_ch = Channel
				.fromPath(params.input_sheet)
        		.splitCsv(sep: '\t', header: ["speci", "genome_id", "genome", "proteins", "genes", "gff", "emapper"])
				.map { gdata -> [gdata.speci, gdata.genome_id, gdata] }

		} else {
			genomes_ch = Channel.fromPath("${params.input_dir}/**")
				.filter( ~/.+\.(ffn|fna|fa(s(ta)?)?)(\.gz)?$/ )
				.map { file -> [ file, file.text.replaceAll(/^>.+$/, "").replaceAll(/\n/, "").length() ] }
				.filter { file, seqlen -> seqlen < params.max_contig_length }
				.map { file, seqlen ->
					def gdata = [:]
					gdata.speci = "contig"
					// gdata.genome_id = file.name.replaceAll(/\.(fasta|fna|fa|ffn)(\.gz)?$/, "")
					gdata.genome_id = fasta.name.replaceAll(/\.(ffn|fna|fa(s(ta)?)?)(\.gz)?$/, "")
					gdata.genome = file
					return [ gdata.speci, gdata.genome_id, gdata ]
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
        		.splitCsv(sep: '\t', header: ["speci", "genome_id", "genome", "proteins", "genes", "gff", "emapper"])
				.map { gdata -> [gdata.speci, gdata.genome_id, gdata] }
        
			speci_ch = genomes_ch
				.map { speci, genome_id, gdata -> speci }
				.unique()
						
		} else {
			// Input genomes are genomic fasta files (.fa, .fasta, .fna, with or without .gz) in a directory or directory tree
			// genomes_ch emits tuples (specI, genome_id, genome_fasta)	
			genomes_ch = Channel.fromPath("${params.input_dir}/**")
				.filter( ~/.+\.(fna|fa(sta)?)(\.gz)?$/ )
				.map { fasta -> 
					def gdata = [:]
					gdata.speci = speci_tag
					gdata.genome_id = fasta.name.replaceAll(/\.(fna|fa(sta)?)(\.gz)?$/, "")
					gdata.genome = fasta
					return [ gdata.speci, gdata.genome_id, gdata ]
				}

			speci_ch = Channel.of(speci_tag)
		}

		status_ch = genomes_ch.map { speci, genome_id, gdata ->
			def flags = [
				GENOME_ANNOTATION: (gdata != null && gdata.genes != null && gdata.proteins != null && gdata.gff != null),
				SPECIES_RECOGNITION: (speci != null && speci != "unknown"),
				SPECI_CLUSTER_SEQS: false,
				RECOMBINASE_SCAN: false,
				FUNCTIONAL_ANNOTATION: (gdata != null && gdata.emapper != null),
				SECRETION_ANNOTATION: false,
				PANGENOME_ESTIMATION: false,
				MGE_ANNOTATION: false
			]
			return [ speci, genome_id, flags ]
		}


		genomes_ch
			.branch {
				speci_annotated: it[0] != "unknown" && it[2].genes != null && it[2].proteins != null && it[2].gff != null // precomputed gene annotation of known species -> recombinase_scan
				speci_unannotated: it[0] != "unknown"                       // genome of known species -> gene_annotation(prodigal)
				speci_unknown: true											// will be determined in species_recognition				
			}
			.set { genomes_speci_ch }

	emit:
		to_recombinase_scan = genomes_speci_ch.speci_annotated
		to_genome_annotation = genomes_speci_ch.speci_unannotated
		to_species_recognition = genomes_speci_ch.speci_unknown
			.map { speci, genome_id, gdata -> [ genome_id, gdata ] }
		speci = speci_ch
		status = status_ch


}
