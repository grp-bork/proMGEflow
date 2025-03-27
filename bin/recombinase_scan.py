#!/usr/bin/env python
# pylint: disable=R0914,E0401


""" Module for parsing recombinase hmmscan results. """

import argparse
import csv
import os

import pyhmmer

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



RECOMBINASE_SCAN_HEADER = (
    "#unigene",
    "recombinase_SMART_hmm_name",
    "PFAM_accession",
    "MGE_prediction",
    "hmmsearch_fullsequence_evalue",
    "hmmsearch_fullsequence_score",
    "MGE_prediction_confidence",
)


def main():
    """ main """
    ap = argparse.ArgumentParser()
    ap.add_argument("proteins", type=str)
    ap.add_argument("hmms", type=str)
    ap.add_argument("--mge_rules", type=str, default=None)
    ap.add_argument("--prefix", type=str, default="sample")
    ap.add_argument("--cpus", "-t", type=int, default=1)
    args = ap.parse_args()

    with pyhmmer.easel.SequenceFile(args.proteins, digital=True) as seq_file:
        protein_seqs = list(seq_file)
    with pyhmmer.plan7.HMMFile(args.hmms) as hmm_file:
        hmm_hits = list(
            pyhmmer.hmmsearch(hmm_file, protein_seqs, bit_cutoffs="gathering")
        )

    raw_table_out = open(
        f"{args.prefix}.recombinase_hmmsearch.out",
        "wb"
    )
    filtered_table_out = open(
        f"{args.prefix}.recombinase_hmmsearch.besthits.out",
        "wb"
    )

    with raw_table_out, filtered_table_out:
        seen = {}
        for i, hits in enumerate(hmm_hits):
            write_header = i == 0
            hits.write(raw_table_out, header=write_header)
            for hit in hits:
                hit_name = hit.name.decode()
                for domain in hit.domains:
                    best_score, _, _ = seen.setdefault(hit_name, (0.0, None, None,))
                    if hit.score > best_score:
                        seen[hit_name] = hit.score, domain, hit
            hits.write(filtered_table_out, header=write_header)

    if seen and args.mge_rules and os.path.isfile(args.mge_rules):
        mge_rules = read_mge_rules(args.mge_rules, recombinase_scan=True)

        with open(
            f"{args.prefix}.recombinase_scan.tsv",
            "wt",
            encoding="UTF-8",
        ) as rscan_out:
            print(*RECOMBINASE_SCAN_HEADER, sep="\t", file=rscan_out)

            for protein, (score, domain, hit) in seen.items():
                hmm_name = domain.alignment.hmm_name.decode()
                print(protein, score, hmm_name)

                recombinase = hmm_name.lower()
                for name, alias in MGE_ALIASES.items():
                    recombinase = recombinase.replace(name, alias)

                rule = mge_rules.get(recombinase)
                if not rule:
                    raise ValueError(f"Cannot find rule for {recombinase}.")

                mges = rule.get_signals()
                confidence = ("low", "high")[len(mges) == 1]

                print(
                    protein,
                    recombinase,
                    domain.alignment.hmm_accession.decode(),
                    ";".join(mges),
                    hit.evalue,
                    hit.score,
                    confidence,
                    sep="\t",
                    file=rscan_out,
                )


if __name__ == "__main__":
    main()
