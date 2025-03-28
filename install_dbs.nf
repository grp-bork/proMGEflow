params.db_dir = "./promge_db"

process setup_emapper_db {
	publishDir "${params.db_dir}", mode: "copy"

	output:
	// tuple val("emapper-5.0.2"), path("eggnog.db"), path("eggnog_proteins.dmnd"), path("eggnog.taxa*"), emit: eggnog_db
	path("emapper-5.0.2"), emit: eggnog_db

	script:
	"""
	set -e

	mkdir -p emapper-5.0.2

	wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
	wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
	wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
	gunzip eggnog.db.gz & 
	gunzip eggnog_proteins.dmnd.gz & 
	tar xvzf eggnog.taxa.tar.gz &
	wait
	rm -vf eggnog.taxa.tar.gz

	mv -v eggnog.db eggnog_proteins.dmnd eggnog.taxa* emapper-5.0.2/
	"""
}

process setup_txsscan_models {
	publishDir "${params.db_dir}", mode: "copy"

	output:
	path("txsscan_models"), emit: txsscan_models

	script:
	"""
	mkdir -p txsscan_models && cd txsscan_models
	git clone https://github.com/macsy-models/TXSScan.git TXSS || \
	(wget https://github.com/macsy-models/TXSScan/archive/refs/heads/master.zip && unzip master.zip && mv -v TXSScan-master TXSS)
	"""
}

process write_param_file {

	publishDir "${params.db_dir}", mode: "copy"

	input:
	path(db_dirs)

	output:
	path("params.yml")

	script:
	def abspath_dbdir = file(params.db_dir).toAbsolutePath()
	"""
	echo -e emapper: >> params.yml
	echo -e "  # path to emapper database" >> params.yml
	echo -e "  db: ${abspath_dbdir}/\$(ls -d emapper*)" >> params.yml
	echo -e txsscan: >> params.yml
	echo -e "  # path to macsyfinder models" >> params.yml
	echo -e "  db: ${abspath_dbdir}/txsscan_models" >> params.yml
	"""

}


workflow install_dbs {

	setup_emapper_db()
	setup_txsscan_models()

	write_param_file(
		setup_emapper_db.out.eggnog_db
			.mix(setup_txsscan_models.out.txsscan_models)
			.collect()
	)
	
}

workflow {
	install_dbs()
}