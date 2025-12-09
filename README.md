# proMGEflow: recombinase-based detection of mobile genetic elements in prokaryotes

`proMGEflow` can detect and annotate mobile genetic elements (MGEs) in prokaryotic genomes based on the presence of recombinase signatures.

![proMGE_workflow](docs/img/proMGEflow.svg)

Dependencies
------------


`proMGEflow` is implemented in `nextflow` and `python`. Dependencies are available as docker containers

- nextflow
- python
- prodigal
- mmseqs2
- eggnog-mapper
- hmmer
- MACsyfinder
- MGExpose
- reCOGnise

Installation
------------

### Databases

`proMGEflow` requires the following databases:

1. eggnog-mapper database (48GB)
	```
	mkdir -p /path/to/emapper_db && cd /path/to/emapper_db
	wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
	wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
	wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
	gunzip eggnog.db.gz
	gunzip eggnog_proteins.dmnd.gz
	tar xvzf eggnog.taxa.tar.gz
	```

2. conjscan models (50MB)
	```
	mkdir -p /path/to/conjscan_models && cd /path/to/conjscan_models
	git clone https://github.com/macsy-models/CONJScan.git CONJ
	cd CONJ
	git checkout d5fc1e3724362cb14c03a6e2f6de879bbdf3f64e
	# ignore the "detached" head message
	```

3. recombinase HMMs (1.7MB)
	```
	wget https://zenodo.org/records/15829523/files/promge_v1_recombinase_models.hmm.gz
	```

4. recognise marker set -- from Zenodo
	```
	TBD
	```

5. pangenome cluster reference sequences (38.9GB)
	```
	wget https://zenodo.org/records/17704403/files/sp095_refdb_v1ypg3.tar
	# Do not untar the .tar file!
	```

Usage
-----

```
nextflow run grp-bork/promgeflow --input_dir /path/to/input/genome/fasta/files --output_dir /path/to/output
```

### Input

#### Via input directory

`proMGEflow` takes as input a set of input genome fasta files. The input files must be stored / linked into a directory, which can be supplied to the workflow via the `--input_dir` parameter. Input genome files can be gzipped and should have a common file ending (supported are `.fa`, `.fasta`, and `.fna`, with or without `.gz` suffix).

#### Via samplesheet

Alternatively, you can supply a samplesheet via the `--input_sheet` parameter. This must be a tab-separated text file with the following columns (without headers):

1. specI cluster id or "unknown" if not known
2. genome id
3. path to genome fasta (.fna)

If gene annotations are already available for the input genome, you can supply them via columns 4-6:

4. path to protein fasta (.faa)
5. path to gene fasta (.ffn)
6. path to (prodigal) gff (.gff)

Additionally, precomputed functional annotations can be specified as well:

7. eggnog-mapper output (.tsv)

### Output

Upon successful MGE detection in an input genome `proMGEflow` returns a gff with the annotated MGEs, the extracted sequences of the MGEs as well as a set of gene coordinates, gene and protein sequences as predicted by `prodigal`.

In case the input genome could not undergo pangenome analysis or no MGEs could be found, `proMGEflow` will return a set of predicted recombinases in the input genome.

