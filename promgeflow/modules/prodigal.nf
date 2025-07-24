process prodigal {
	tag "${genome_id}"
	label "prodigal"
	container "quay.io/biocontainers/prodigal:2.6.3--h031d066_7"
	cpus 1
	time {1.d * task.attempt}
	memory {8.GB * task.attempt}

	input:
	tuple val(speci), val(genome_id), path(genome_fna)

	output:
	tuple val(speci), val(genome_id), path("${genome_id}/${genome_id}.faa"), emit: proteins
	tuple val(speci), val(genome_id), path("${genome_id}/${genome_id}.ffn"), emit: genes
	tuple val(speci), val(genome_id), path("${genome_id}/${genome_id}.gff"), emit: genome_annotation

	script:
	def gunzip_cmd = (genome_fna.name.endsWith(".gz")) ? "gzip -dc ${genome_fna} > \$(basename ${genome_fna} .gz)" : ""
	"""
	mkdir -p ${genome_id}
	${gunzip_cmd}
	prodigal -i \$(basename ${genome_fna} .gz) -f gff -o ${genome_id}/${genome_id}.gff -a ${genome_id}/${genome_id}.faa -d ${genome_id}/${genome_id}.ffn
	"""
}

process buffered_prodigal {
	tag "batch_${batch_id}"
	label "prodigal"
	container "quay.io/biocontainers/prodigal:2.6.3--h031d066_7"
	cpus 1
	time {1.d * task.attempt}
	memory {8.GB * task.attempt}

	input:
	tuple val(batch_id), path(genomes)
	val(file_suffix)

	output:
	path("prodigal/**.*"), emit: annotations

	script:
	"""
	for genome_file in ${genomes}; do
		# genome_id=\$(basename \$genome_file ${file_suffix})
		genome_id=\$(basename \$genome_file)
		echo \$genome_id
		# mkdir -p prodigal/\$genome_id/
		if [[ \$genome_file == *.gz ]]; then
			gzip -dc \$genome_file > \$(basename \$genome_file .gz) 			
		fi
		# prodigal -i \$(basename \$genome_file .gz) -f gff -o prodigal/\$genome_id/\$genome_id.gff -a prodigal/\$genome_id/\$genome_id.faa -d prodigal/\$genome_id/\$genome_id.ffn
		prodigal -i \$(basename \$genome_file .gz) -f gff -o prodigal/\$genome_id.gff -a prodigal/\$genome_id.faa -d prodigal/\$genome_id.ffn
		if [[ \$genome_file == *.gz ]]; then
			rm -fv \$(basename \$genome_file .gz) 			
		fi
	done
	"""

}

process publish_annotations {
	publishDir "${params.output_dir}", mode: "copy"
	executor "local"
	tag "${genome_id}"

	input:
	tuple val(speci), val(genome_id), path(annotations)

	output:
	path("${speci}/${genome_id}/**")

	script:
	"""
	mkdir -p ${speci}/${genome_id}/

	ln -s ../../${genome_id}.faa ${speci}/${genome_id}/
	ln -s ../../${genome_id}.ffn ${speci}/${genome_id}/
	ln -s ../../${genome_id}.gff ${speci}/${genome_id}/
	"""
}
