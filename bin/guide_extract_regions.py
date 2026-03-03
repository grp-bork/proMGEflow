#!/usr/bin/env python3

import csv
import gzip
import sys




def get_lines_from_chunks(f: str, bufsize: int = 800000000):
    """
    Provides generator access to the lines of large text files.
    File is read chunk-wise into a buffer of the specified size.
    Support gzip-compressed files.

    inputs:
     - f: str -- filename
     - bufsize: int -- size of buffer

    """
    gzipped = f.endswith(".gz")
    with (gzip.open if gzipped else open)(f, "r") as _in:
        tail = ""
        while 1:
            chunk = _in.read(bufsize)
            if gzipped:
                chunk = chunk.decode()
            chunk = "".join((tail, chunk))
            if not chunk:
                break
            chunk = chunk.split("\n")
            *lines, tail = chunk
            yield from lines
        if tail:
            yield tail


def read_fasta(f):
    header, seq = None, []
    for line in get_lines_from_chunks(f):
        if line[0] == ">":
            if seq:
                yield header, "".join(seq)
                seq.clear()
            header = line.strip()[1:]
        else:
            seq.append(line.strip())
    if seq:
        yield header, "".join(seq)


def process_island(island, genes, genome_id, stream=sys.stdout,):

    # island_id_match = re.match(r'(GIL|MGE|SPIRE)_(.+)_(.+):(\d+)-(\d+)', id_string)
    #     if not island_id_match:
    #         raise ValueError(f"{id_string} does not seem to be a valid island identifier.")

    #     groups = island_id_match.groups()[1:]
    
    start = min(int(gene[3]) for gene in genes)
    end = max(int(gene[4]) for gene in genes)
    island_id = f"GIL_{genome_id}_{island[0]}:{start}-{end}"

    print(
        island[0],
        "proMGE_guide",
        "region",
        # int(island[1]) + 1,
        start,
        # island[2],
        end,
        ".",
        ".",
        ".",
        f"ID={island_id};{island[3]}",
        file=stream,
        sep="\t",
    )
    for gene in sorted(genes, key=lambda x:tuple(map(int, x.split("\t")[3:5]))):
        gene[2] = "gene"
        print(*gene, sep="\t", file=stream,)


def main():
    # k119_100010     266     9537    recombinase=xer,355-1281;n_aln=4;hc_region=269-9163,0.785;lc_cov=0.763  k119_100010     Prodigal_v2.6.3 CDS     355     1281    57.5    -       0       ID=143600_1;partial=00;start_type=ATG;rbs_motif=AACAA;rbs_spacer=14bp;gc_cont=0.384;con>
    # k119_100010     Prodigal_v2.6.3 CDS     355     1281    57.5    -       0       ID=143600_1;partial=00;start_type=ATG;rbs_motif=AACAA;rbs_spacer=14bp;gc_cont=0.384;con>
    island = None
    genes = set()

    gene_coords = {}

    with open(sys.argv[1], 'rt') as _in, open(f"{sys.argv[3]}.mge_candidates.gff3", "wt") as gff_out:
        print("##gff-version 3", file=gff_out,)

        for row in csv.reader(_in, delimiter='\t'):
            new_island = row[:3]  # was :4 but then we have duplicate islands in case of tandem recombinases
            if island != new_island:
                if island is not None:
                    process_island(island, genes, sys.argv[3], stream=gff_out,)
                    # print(
                    #     island[0],
                    #     "proMGE_guide",
                    #     "region",
                    #     int(island[1]) + 1,
                    #     island[2],
                    #     ".",
                    #     ".",
                    #     ".",
                    #     island[3],
                    #     file=gff_out,
                    #     sep="\t",
                    # )
                    # for gene in genes:
                    #     gene[2] = "gene"
                    #     print(*gene, sep="\t", file=gff_out,)

                island = new_island
                genes.clear()

            genes.add(row[4:-1])
            gene_coords.setdefault(genes[-1][0], set()).add(tuple(map(int, genes[-1][3:5])))

        if island is not None:
            process_island(island, genes, sys.argv[3], stream=gff_out,) 
            # print(
            #     island[0],
            #     "proMGE_guide",
            #     "region",
            #     int(island[1]) + 1,
            #     island[2],
            #     ".",
            #     ".",
            #     ".",
            #     island[3],
            #     file=gff_out,
            #     sep="\t",
            # )
            # for gene in genes:
            #     gene[2] = "gene"
            #     print(*gene, sep="\t", file=gff_out,)

    with open(f"{sys.argv[3]}.cargo.faa", "wt") as faa_out:
        for sid, seq in read_fasta(sys.argv[2]):
            # >k119_155575_1 # 3 # 230 # -1 # ID=1_1;partial=10;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.509
            tokens = sid.strip().split(" # ")
            contig = tokens[0][:tokens[0].rfind("_")]
            start, end = map(int, tokens[1:3])

            if (start, end) in gene_coords.get(contig, set()):
                print(f">{sid}\n{seq}", file=faa_out,)





if __name__ == "__main__":
    main()