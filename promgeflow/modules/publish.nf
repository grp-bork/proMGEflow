process publish_results {
	label "tiny"
	tag "${speci}/${genome_id}"

	input:
	tuple val(speci), val(genome_id), path("promgeflow_results/*")
	val(simple_output)
	val(as_tarball)

	output:
	path("**.*")

	script:
	def outdir = "${speci}/${genome_id}"
	def lvlup = "../.."

	if (as_tarball) {
		"""
		tar cvzf ${genome_id}.promgeflow.tar.gz promgeflow_results/ 
		"""
	} else {
		if (simple_output) {
			outdir = "${genome_id}"
			lvlup = ".."
		}
		"""
		mkdir ${outdir} && cd ${outdir}
		find ${lvlup} -type l -exec ln -s {} \\;
		"""
	}
		
}


process publish_gene_annotations {
	// executor "local"  -> move to run.config @ EMBL
	label "tiny"
	tag "${speci}/${genome_id}"

	input:
	tuple val(speci), val(genome_id), path(annotations)
	val(simple_output)

	output:
	path("**/*.{faa,ffn,gff}")

	script:

	def outdir = "${speci}/${genome_id}"
	def lvlup = "../.."

	if (simple_output) {
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
	// executor "local"  -> move to run.config @ EMBL
	label "tiny"
	tag "${speci}/${genome_id}"

	input:
	tuple val(speci), val(genome_id), path(mge_table), path(mge_gff)
	val(simple_output)

	output:
	path("**/*.{tsv,gff3}")

	script:
	def outdir = "${speci}/${genome_id}"
	def lvlup = "../.."
	if (simple_output) {
		outdir = "${genome_id}"		
		lvlup = ".."
	}

	"""
	mkdir -p ${outdir}/ && cd ${outdir}

	ln -s ${lvlup}/*.tsv
	ln -s ${lvlup}/*.gff3
	"""
}