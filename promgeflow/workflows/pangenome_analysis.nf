include { linclust } from "../modules/linclust"


workflow pangenome_analysis {

	take:
		genomes_ch
		cluster_reps_ch

	main:
		genomes_ch.dump(pretty: true, tag: "pangenome_analysis_input")

		genes_ch = genomes_ch
			.filter { it[3].FUNCTIONAL_ANNOTATION }
			.map { speci, genome_id, gdata, flags -> [ speci, genome_id, gdata.genes ] }

		linclust_input_ch = genes_ch
			.combine(cluster_reps_ch, by: 0)

		linclust(linclust_input_ch)

		linclust_clusters_ch = linclust.out.mmseq_cluster
			.join(linclust.out.done_sentinel, by: [0, 1], remainder: true)
			.map { speci, genome_id, clusters, sentinel ->
				[ speci, genome_id, (sentinel != null) ? clusters : null ]
			}
		linclust_clusters_ch = genomes_ch
			.join(linclust_clusters_ch, by: [0, 1], remainder: true)
			.map { speci, genome_id, gdata_old, flags_old, clusters ->
				def gdata = gdata_old.clone()
				gdata.gene_clusters = clusters
				def flags = flags_old.clone()
				flags.PANGENOME_CLUSTERING = (clusters != null)
				return [ speci, genome_id, gdata, flags ]
			}

		linclust_clusters_ch.dump(pretty: true, tag: "linclust_clusters_ch")
	
	emit:
		genomes = linclust_clusters_ch

}