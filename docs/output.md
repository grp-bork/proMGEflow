# Output

`proMGEflow` produces the following outputs:

## Directory Structure and Main Output Files

Genomes, for which a specI was supplied (or identified) will be in the following output directory tree:

- `<specI>/<genome>`  

  - `<genome>.faa` - prodigal protein sequences
  - `<genome>.ffn` - prodigal gene sequences
  - `<genome>.gff` - prodigal gene predictions
  - `<genome>.mge_islands.gff3` - proMGE/mgexpose MGE annotations
  - `<genome>.mge_islands.ffn.gz` - DNA sequences of the annotated MGEs

Genomes without assigned specI will be stored in

- `unknown/<genome>`

  - `<genome>.MGE_predictions.tsv` - table with hits to proMGEs recombinase HMM database