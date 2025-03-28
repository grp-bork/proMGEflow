include { full_annotation } from "./promgeflow/workflows/full_annotation"
include { install_dbs } from "./install_dbs"

workflow {

	if (params.install_dbs) {
		install_dbs()
	} else {
		full_annotation()
	}



}