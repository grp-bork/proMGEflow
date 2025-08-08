#!/usr/bin/env python3

""" Module for parsing recombinase hmmscan results. """


import argparse
import csv
import re

from recombinases import MGE_ALIASES, MgeRule

def read_mge_rules(f, recombinase_scan=False):
    """ Read MGE rules.

    Returns dictionary {mge: MgeRule}.
    """
    with open(f, "rt", encoding="UTF-8") as _in:
        rules = {
            row[0].lower(): MgeRule(row[0], *(tuple(map(int, row[1:]))), recombinase_scan)
            for i, row in enumerate(csv.reader(_in, delimiter="\t"))
            if i != 0
        }

    # #special case for Tn3 since it can carry conjugative system#
    # for rule_id, rule in rules.items():
    # 	if "tn3" in rule_id:
    # 		rule.ce = 1

    return rules


def parse_hmm_table(table_stream):
    """ Parse hmm table. """
    for line in table_stream:
        line = line.strip()
        if line and line[0] != "#":
            parsed = re.split(r"\s+", line)
            yield parsed, parsed[0], float(parsed[5])


def extract_best_hits(table_stream):
    """ Extracts best hits from hmm-table stream.

    Returns best-scoring recombinase hmm hits.
    """
    seen = {}
    for line, protein, score in table_stream:
        seen_score = seen.setdefault(protein, [0.0, ""])[0]
        if score > seen_score:
            seen[protein] = [score, line]

    return [line for _, line in seen.values()]


def generate_output_table(best_hits, mge_rules):
    """ Annotates recombinase hmm hits using mge rules.

    Returns annotated table rows via generator.
    """
    for hit in best_hits:
        hit[2] = hit[2].lower()
        for name, alias in MGE_ALIASES.items():
            hit[2] = hit[2].replace(name, alias)

        rule = mge_rules.get(hit[2])
        if not rule:
            raise ValueError(f"Cannot find rule for {hit[2]}.")

        mges = rule.get_signals()
        confidence = "high" if len(mges) == 1 else "low"

        yield (hit[0], hit[2], hit[3], ";".join(mges), hit[4], hit[5], confidence)


def main():
    """ Main function, duh. """
    ap = argparse.ArgumentParser()
    ap.add_argument("hmmsearch_table", type=str)
    ap.add_argument("--mge_rules", type=str, default=None)
    ap.add_argument("--prefix", type=str, default="sample")
    ap.add_argument("--proteins", type=str)
    args = ap.parse_args()

    proteins = {}
    if args.proteins:
        with open(args.proteins, "rt") as _in:
            for line in _in:
                if line[0] == ">":
                    #  >gnl|AGSH|NT12270_11_4 # 2457 # 3029 # 1 # ID=11_4;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.483
                    line = line.strip().split(" ")
                    proteins[line[0][1:]] = {
                        "start": line[2],
                        "end": line[4],
                        "strand": "+" if line[6] == "1" else "-",
                        "attribs": line[8],
                    }

    best_hits = []
    with open(args.hmmsearch_table, "rt", encoding="UTF-8") as table_stream:
        best_hits = extract_best_hits(parse_hmm_table(table_stream))

    if best_hits:
        with open(
            f"{args.prefix}.recombinase_hmmsearch.besthits.out",
            "wt",
            encoding="UTF-8",
        ) as raw_table_out:
            print(*("\t".join(bh) for bh in best_hits), sep="\n", file=raw_table_out)

        if args.mge_rules:
            mge_rules = read_mge_rules(args.mge_rules, recombinase_scan=True)

            mge_pred_out = open(
                f"{args.prefix}.recombinase_based_MGE_predictions.tsv",
                "wt",
                encoding="UTF-8",
            )

            mge_pred_gff = open(
                f"{args.prefix}.predicted_recombinase_mges.gff3",
                "wt",
                encoding="UTF-8",
            )

            with mge_pred_out, mge_pred_gff:

                header = (
                    "#unigene", "recombinase_SMART_hmm_name", "PFAM_accession",
                    "MGE_prediction", "hmmsearch_fullsequence_evalue",
                    "hmmsearch_fullsequence_score", "MGE_prediction_confidence"
                )

                print(*header, sep="\t", file=mge_pred_out)
                print("##gff-version 3", file=mge_pred_gff)

                recombinases = []
                for line in generate_output_table(best_hits, mge_rules):
                    print(*line, sep="\t", file=mge_pred_out)
                    protein = proteins.get(line[0])
                    if protein is not None:
                        attribs=";".join(f"{k}={v.replace(';', ',')}" for k, v in zip(("recombinase", "PFAM", "predicted_mge", "evalue", "score", "confidence",), line[1:]))
                        
                        rline = (
                            line[0][:line[0].rfind("_")],
                            "proMGE",
                            "gene",
                            protein["start"],
                            protein["end"],
                            protein["strand"],
                            ".",
                            ";".join((attribs, protein["attribs"],))
                        )
                        recombinases.append(rline)

                for line in sorted(recombinases, key=lambda x:(x[0], int(x[3]), int(x[4]))):
                    # gnl|AGSH|NT12270_27_3   dde_tnp_is1     PF03400.12      is_tn   3.1e-74 245.6   high
                    print(*line, sep="\t", file=mge_pred_gff,)

                        


if __name__ == "__main__":
    main()
