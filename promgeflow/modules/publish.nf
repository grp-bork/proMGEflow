process publish_tarball {
	label "tiny"
	// tag "${speci}/${genome_id}"
	tag "Publishing results..."

	input:
	// tuple val(speci), val(genome_id), path("promgeflow_results/*")
	path("promgeflow_results_raw/*")
	// val(simple_output)
	val(as_tarball)

	output:
	// path("**.*")
	path("*.tar.gz")

	script:
	// def outdir = "${speci}/${genome_id}"
	def lvlup = "../.."

	if (true) {
		def tarball_prefix = (as_tarball && as_tarball?.trim()) ? "${as_tarball}" : "promgeflow"
		"""
		mkdir -p promgeflow_results/

		for f in \$(find promgeflow_results_raw -name '*.gff3'); do
			s=\$(basename \$f | sed "s/\\.\\(mge_islands\\|predicted_recombinase_mges\\)\\.gff3//");
			mkdir -p promgeflow_results/\$s;
			find promgeflow_results_raw -name '\$s*' -exec ln -sf ../../{} promgeflow_results/\$s \\;
		done

		tar chvzf ${tarball_prefix}.tar.gz promgeflow_results/ 
		"""
		//  | xargs -I sh -c 's=\$(basename {} | sed "s/\\.\\(mge_islands\\|predicted_recombinase_mges\\)\\.gff3//""); mkdir -p promgeflow_results/\$s; ln -sf ../../promgeflow_results_raw/\$s* promgeflow_results/\$s/'
	} 
	// else {
	// 	if (simple_output) {
	// 		outdir = "${genome_id}"
	// 		lvlup = ".."
	// 	}
	// 	"""
	// 	mkdir ${outdir} && cd ${outdir}
	// 	find ${lvlup} -type l -exec ln -s {} \\;
	// 	"""
	// }
		
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