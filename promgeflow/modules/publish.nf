process publish_tarball {
	label "tiny"
	tag "Publishing results..."

	input:
	path("promgeflow_results_raw/*")
	val(as_tarball)

	output:
	path("*.tar.gz")

	script:
	def tarball_prefix = (as_tarball && as_tarball?.trim()) ? "${as_tarball}" : "promgeflow_results"
	"""
	mkdir -p ${tarball_prefix}/

	for f in \$(find promgeflow_results_raw -name '*.gff3'); do
		s=\$(basename \$f | sed "s/\\.\\(mge_islands\\|predicted_recombinase_mges\\)\\.gff3//");
		mkdir -p ${tarball_prefix}/\$s;
		find promgeflow_results_raw -name "\$s*" -exec ln -sf ../../{} ${tarball_prefix}/\$s \\;
	done

	find promgeflow_results_raw -name '*.txt' -exec ln -sf ../{} ${tarball_prefix}/ \\;

	find ${tarball_prefix} -name '*.fna.???' | xargs -I{} sh -c 't=\$(ls {} | sed "s/\\.fna\\.\\(faa\\|ffn\\|gff\\)/.\\1/"); mv -v {} \$t;'

	tar chvzf ${tarball_prefix}.tar.gz ${tarball_prefix}/ 
	"""
		
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