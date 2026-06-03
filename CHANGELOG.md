v3.13.5
-------

* update reCOGnise to v0.8

v3.13.4 (no clowm)
------------------

* add missing column to recombinase_scan gff

v3.13.3 (no clowm)
------------------

* fix MGExpose container URL

v3.13.2 (no clowm)
------------------

* fix issue with cluster parsing

v3.13.1 (no clowm)
------------------

* fixed issue in linclust input that would cause only query genes to be clustered
* fixed issue in mgexpose that would lead to false gene and pangenome counts

v3.13.0 (no clowm)
------------------

* updated samplesheet input to support precomputed genome annotations as well as functional gene annotations

* secretion annotation now uses CONJScan models instead of TXSScan models
  * txsscan process was renamed to macsyfinder

* refactored and decluttered channel routing in de novo workflow

* added --force_secretion_analysis parameter (e.g. for bulk annotation)

* added gene_info.txt to mgexpose published outputs

* re-enabled gzip inputs for most processes to support precomputed, gzipped inputs via sample sheet

* re-added and updated gene-level recognise

* added output sentinel to recombinase scan process

* added pangenome reporting