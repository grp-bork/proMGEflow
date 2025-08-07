// If the specI classification is known beforehand, it can be passed via the `known_speci` parameter.
// Otherwise, specI will be treated as unknown until it is detected by reCOGnise downstream.
def speci_tag = params.known_speci ?: "unknown"


workflow handle_input_plasmids {

	main:
		def reg_ctr = 0
    	plasmids_ch = Channel
			.fromPath(params.input_fasta)
			.splitFasta(by: 1, file: true)
			.map { file -> [ file, file.text.replaceAll(/^>.+$/, "").replaceAll(/\n/, "").length() ] }
			.filter { file, seqlen -> seqlen < params.max_plasmid_length }
			.map { file, seqlen ->
				def genome = file.name.replaceAll(/\.[0-9]+\.(fasta|fna|fa|ffn)(\.[2a-z]+)?$/, "")
				def contig = file.text.split("\n")[0].split(" ")[0].substring(1) //replaceAll(/^>([^ ]+).+/, "\1")
				def region_id = "${reg_ctr++}_CP_100.${contig}.1-${seqlen}.${contig}"
				[ genome, contig, region_id, file ]
			}
			// .map { file, seqlen ->
			// 	def contig = file.name.replaceAll(/\.[0-9]+\.(fasta|fna|fa|ffn)(\.[2a-z]+)?$/, "")
			// 	def region_id = "${reg_ctr++}_CP_100.${contig}.1-${seqlen}.${contig}"
			// 	[ contig, region_id, file ]
			// }

		plasmids_ch.dump(pretty: true, tag: "plasmids_ch")

		genomes_ch = plasmids_ch.map { genome_id, contig_id, region_id, file -> [ "plasmid", genome_id, file ] }
		regions_ch = plasmids_ch.map { genome_id, contig_id, region_id, file -> [ "plasmid", genome_id, region_id ] }

	emit:
		genomes = genomes_ch
		regions = regions_ch

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
			// Input genomes are genomic fasta files (.fa, .fasta, .fna, with or without .gz) in a directory or directory tree
			// genomes_ch emits tuples (specI, genome_id, genome_fasta)	
			genomes_ch = Channel.fromPath("${params.input_dir}/**")
				.filter( ~/.+\.(fna|fa(sta)?)(\.gz)?$/ )
				.map { fasta -> [ speci_tag, file.name.replaceAll(/\.(fna|fa(sta)?)(\.gz)?$/, ""), fasta ] }

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