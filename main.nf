include { full_annotation } from "./promgeflow/workflows/full_annotation"
include { contig_annotation } from "./promgeflow/workflows/contig_annotation"


workflow {

	if (params.run_mode == "contig") {
		contig_annotation()
	} else {
		full_annotation()
	}

}