process publish_gene_annotations {
	executor "local"
	tag "${speci}/${genome_id}"

	input:
	tuple val(speci), val(genome_id), path(annotations)
	val(simple_outdir)

	output:
	path("**/*.{faa,ffn,gff}")

	script:

	def outdir = "${speci}/${genome_id}"
	def lvlup = "../.."

	if (simple_outdir) {
		outdir = "${genome_id}"
	 	lvlup = ".."
	}

	"""
	mkdir -p ${outdir}/ && cd ${outdir}

	ln -s ${lvlup}/*.faa ${genome_id}.faa
	ln -s ${lvlup}/*.ffn ${genome_id}.ffn
	ln -s ${lvlup}/*.gff ${genome_id}.gff
	"""
}

process publish_recombinase_scan {
	executor "local"
	tag "${speci}/${genome_id}"

	input:
	tuple val(speci), val(genome_id), path(mge_table), path(mge_gff)
	val(simple_outdir)

	output:
	path("**/*.{tsv,gff3}")

	script:
	def outdir = "${speci}/${genome_id}"
	def lvlup = "../.."
	if (simple_outdir) {
		outdir = "${genome_id}"		
		lvlup = ".."
	}

	"""
	mkdir -p ${outdir}/ && cd ${outdir}

	ln -s ${lvlup}/*.tsv
	ln -s ${lvlup}/*.gff3
	"""
}