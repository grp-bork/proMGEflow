# proMGEflow: recombinase-based detection of mobile genetic elements in prokaryotes

Blahblah.

![proMGE_workflow](https://raw.githubusercontent.com/grp-bork/proMGEflow/main/docs/img/proMGEflow.svg)

Dependencies
------------


`proMGEflow` is implemented in `nextflow` and `python`. Dependencies are available as docker containers

- nextflow
- python
- prodigal
- mmseqs2
- eggnog-mapper
- hmmer
- MACsyfinder/txsscan
- MGExpose
- reCOGnise

Installation
------------

### Databases

`proMGEflow` requires the following databases:

1. eggnog-mapper database --
2. txsscan database --
3. recombinase HMMs -- from Zenodo
4. recognise marker set -- from Zenodo
5. pangenome cluster reference sequences -- from Zenodo

Usage
-----

```
nextflow run grp-bork/promgeflow --input_dir /path/to/input/genome/fasta/files --output_dir /path/to/output
```

### Input

`proMGEflow` takes as input a set of input genome fasta files. The input files must be stored / linked into a directory, which can be supplied to the workflow via the `--input_dir` parameter. Input genome files can be gzipped and should have a common file ending (default `.fna`). The file ending can be supplied to the workflow via the `--file_pattern` parameter (e.g. `--file_pattern '*.fasta'`).

### Output

TBD

