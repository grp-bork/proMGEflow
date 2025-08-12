include { full_annotation } from "./promgeflow/workflows/full_annotation"
include { install_dbs } from "./install_dbs"
include { plasmid_annotation } from "./promgeflow/workflows/plasmid_annotation"


workflow {

	if (params.run_mode == "install_dbs") {
		install_dbs()
	} else if (params.run_mode == "plasmid") {
		plasmid_annotation()
	} else {
		full_annotation()
	}

}