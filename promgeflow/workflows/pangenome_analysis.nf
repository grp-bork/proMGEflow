include { linclust } from "../modules/linclust"


workflow pangenome_analysis {

	take:
		filtered_genes_ch
		cluster_reps_ch

	main:
		linclust_input_ch = filtered_genes_ch
			.combine(cluster_reps_ch, by: 0)

		linclust(linclust_input_ch)

		linclust_clusters_ch = linclust.out.mmseq_cluster
			.join(linclust.out.done_sentinel, by: [0, 1])
			.map { speci, genome_id, clusters, sentinel ->
				[ speci, genome_id, clusters ]
			}		

		linclust_clusters_ch.dump(pretty: true, tag: "linclust_clusters_ch")
	
	emit:
		clusters = linclust_clusters_ch

}