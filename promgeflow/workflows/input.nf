// If the specI classification is known beforehand, it can be passed via the `known_speci` parameter.
// Otherwise, specI will be treated as unknown until it is detected by reCOGnise downstream.


workflow handle_input_plasmids {

	main:
		def reg_ctr = 0
    	plasmids_ch = channel
			.fromPath(params.input_fasta)
			.splitFasta(by: 1, file: true)
			.map { file -> [ file, file.text.replaceAll(/^>.+$/, "").replaceAll(/\n/, "").length() ] }
			.filter { _file, seqlen -> seqlen < params.max_plasmid_length }
			.map { file, seqlen ->
				def genome = file.name.replaceAll(/\.[0-9]+\.(fasta|fna|fa|ffn)(\.[2a-z]+)?$/, "")
				def contig = file.text.split("\n")[0].split(" ")[0].substring(1) //replaceAll(/^>([^ ]+).+/, "\1")
				def region_id = "${reg_ctr}_CP_100.${contig}.1-${seqlen}.${contig}"
				reg_ctr += 1
				return [ reg_ctr, genome, contig, region_id, file ]
			}
			.map { _ctr, genome, contig, region_id, file -> [ genome, contig, region_id, file ] }  // double-map is brutal but only solution for now to make v26 compatible..
			// .map { file, seqlen ->
			// 	def contig = file.name.replaceAll(/\.[0-9]+\.(fasta|fna|fa|ffn)(\.[2a-z]+)?$/, "")
			// 	def region_id = "${reg_ctr++}_CP_100.${contig}.1-${seqlen}.${contig}"
			// 	[ contig, region_id, file ]
			// }

		plasmids_ch.dump(pretty: true, tag: "plasmids_ch")

		genomes_ch = plasmids_ch.map { genome_id, _contig_id, _region_id, file -> [ "plasmid", genome_id, file ] }
		regions_ch = plasmids_ch.map { genome_id, _contig_id, region_id, _file -> [ "plasmid", genome_id, region_id ] }

	emit:
		genomes = genomes_ch
		regions = regions_ch

}


workflow handle_input_genomes {		

	main:
		def speci_tag = params.known_speci ?: "unknown"
		genomes_ch = channel.empty()
		speci_ch = channel.empty()

		if (params.input_sheet) {
			genomes_ch = channel
				.fromPath(params.input_sheet)
        		.splitCsv(sep: '\t', header: ["speci", "genome_id", "genome", "proteins", "genes", "gff", "emapper"])
				.map { gdata -> [gdata.speci, gdata.genome_id, gdata] }
        		// .map { it -> [it[0], it[1], it[2]] }
        
			speci_ch = genomes_ch
				.map { speci, _genome_id, _gdata -> speci }
				// .map {speci, genome, file -> speci}
				.unique()
				.view()
						
		} else {
			// Input genomes are genomic fasta files (.fa, .fasta, .fna, with or without .gz) in a directory or directory tree
			// genomes_ch emits tuples (specI, genome_id, genome_fasta)	
			genomes_ch = channel.fromPath("${params.input_dir}/**")
				.filter( ~/.+\.(fna|fa(sta)?)(\.gz)?$/ )
				.map { fasta -> 
					def gdata = [:]
					gdata.speci = speci_tag
					gdata.genome_id = fasta.name.replaceAll(/\.(fna|fa(sta)?)(\.gz)?$/, "")
					gdata.genome = fasta
					return [ gdata.speci, gdata.genome_id, gdata ]
				}

			speci_ch = channel.of(speci_tag)
		}

		genomes_ch
			.branch {
				// speci_known: it[0] != "unknown"
				// speci_unknown: true
				speci_annotated: { speci, _genome_id, gdata -> speci != "unknown" && gdata.genes != null }  // precomputed gene annotation of known species -> recombinase_scan
				speci_unannotated: { speci, _genome_id, _gdata -> speci != "unknown" }                       // genome of known species -> gene_annotation(prodigal)
				speci_unknown: true
				// will be determined in species_recognition
				// unknown_annotated: it[2].genes != null                      // precomputed gene annotation of unknown species -> recognise(_genes)
				// unknown_unannotated: true                                   // genome of unknown species -> recognise_genome
			}
			.set { genomes_speci_ch }


	emit:
		to_recombinase_scan = genomes_speci_ch.speci_annotated
		to_genome_annotation = genomes_speci_ch.speci_unannotated
		to_species_recognition = genomes_speci_ch.speci_unknown
			.map { _speci, genome_id, gdata -> [ genome_id, gdata ] }
		// genomes_with_speci = genomes_speci_ch.speci_known
		// genomes_without_speci = genomes_speci_ch
		// 	.speci_unknown
		// 	.map { speci, genome_id, gdata -> [genome_id, gdata] }
		speci = speci_ch


}