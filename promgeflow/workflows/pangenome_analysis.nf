include { linclust } from "../modules/linclust"


workflow pangenome_analysis {

	take:
		// filtered_genes_ch
		genomes_ch
		cluster_reps_ch

	main:

		genes_ch = genomes_ch.map { speci, genome_id, gdata -> [ speci, genome_id, gdata.genes ] }

		linclust_input_ch = genes_ch
			.combine(cluster_reps_ch, by: 0)

		linclust(linclust_input_ch)

		linclust_clusters_ch = linclust.out.mmseq_cluster
			.join(linclust.out.done_sentinel, by: [0, 1])
			.map { speci, genome_id, clusters, sentinel ->
				[ speci, genome_id, clusters ]
			}
			.join(genomes_ch)
			.map { speci, genome_id, clusters, gdata_old ->
				def gdata = gdata_old.clone()
				gdata.gene_clusters = clusters
				return [ speci, genome_id, clusters, gdata ]
			}	

		linclust_clusters_ch.dump(pretty: true, tag: "linclust_clusters_ch")
	
	emit:
		genomes = linclust_clusters_ch

}