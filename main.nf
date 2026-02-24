include { full_annotation } from "./promgeflow/workflows/full_annotation"
include { plasmid_annotation } from "./promgeflow/workflows/plasmid_annotation"
include { guided_annotation } from "./promgeflow/workflows/guided"


workflow {

	if (params.run_mode == "plasmid") {
		plasmid_annotation()
	} else if (params.run_mode == "guided") {
		guided_annotation()
	} else {
		full_annotation()
	}

}