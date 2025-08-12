// phasing out nested parameters
params.recognise = [:]
params.recognise.marker_set = "motus"
params.recognise_marker_set = params.recognise.marker_set


process recognise_genome {
	container "ghcr.io/grp-bork/recognise:main"
	tag "${genome_id}"
	label "recognise"
	label "small"

	input:
	tuple val(genome_id), path(genome)
	path(marker_genes_db)
	
	
	output:
	tuple val(genome_id), path("**/${genome_id}.specI.status.OK"), emit: speci_status_ok, optional: true 
	tuple val(genome_id), path("**/${genome_id}.cogs.txt"), emit: cog_table
	tuple val(genome_id), path("**/${genome_id}.specI.txt"), emit: genome_speci
	tuple val(genome_id), path("**/${genome_id}.specI.status"), emit: speci_status
	tuple val(genome_id), path("**/${genome_id}.faa"), emit: proteins
	tuple val(genome_id), path("**/${genome_id}.ffn"), emit: genes
	tuple val(genome_id), path("**/${genome_id}/${genome_id}.gff"), emit: gff
	
	script:
	"""
	if [[ "${genome}" == *".gz" ]]; then
		gzip -dc ${genome} > genome_file
	else
		ln -sf ${genome} genome_file
	fi

	recognise --marker_set ${params.recognise_marker_set} --genome genome_file --cpus ${task.cpus} --with_gff -o recognise/${genome_id} ${genome_id} \$(readlink ${marker_genes_db})

	if [[ -s recognise/${genome_id}/${genome_id}.specI.txt ]]; then 
		speci=\$(cat recognise/${genome_id}/${genome_id}.specI.txt)
	else
		speci="unknown"
	fi

	mv -v recognise \$speci
	
	rm -fv genome_file
	"""

}
