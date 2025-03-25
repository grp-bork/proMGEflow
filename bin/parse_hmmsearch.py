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
    args = ap.parse_args()

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

            with open(
                f"{args.prefix}.recombinase_based_MGE_predictions.tsv",
                "wt",
                encoding="UTF-8",
            ) as mge_pred_out:

                header = (
                    "#unigene", "recombinase_SMART_hmm_name", "PFAM_accession",
                    "MGE_prediction", "hmmsearch_fullsequence_evalue",
                    "hmmsearch_fullsequence_score", "MGE_prediction_confidence"
                )

                print(*header, sep="\t", file=mge_pred_out)

                for line in generate_output_table(best_hits, mge_rules):
                    print(*line, sep="\t", file=mge_pred_out)


if __name__ == "__main__":
    main()
