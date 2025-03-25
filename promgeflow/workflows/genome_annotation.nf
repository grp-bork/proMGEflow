include { prodigal } from "../modules/prodigal"


workflow genome_annotation {

	take:
		genomes_ch
	
	main:
	    // suffix_pattern and speci_tag are defined in main.nf but not recognized otherwise
	    def suffix_pattern = params.file_pattern.replaceAll(/\*/, "")
	    def speci_tag = params.known_speci ?: "unknown_speci"

		prodigal(genomes_ch)
		pproteins_ch = prodigal.out.proteins
		pgenes_ch = prodigal.out.genes
		pgffs_ch = prodigal.out.genome_annotation

	emit:
		proteins = pproteins_ch
		genes = pgenes_ch
		gffs = pgffs_ch

}