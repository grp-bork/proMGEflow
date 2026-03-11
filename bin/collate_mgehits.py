#!/usr/bin/env python3

import csv
import glob
import sys

def main():

	d = {}

	for i, f in enumerate(glob.glob(f"{sys.argv[1]}/*.tsv")):
		with open(f, "rt") as _in:
			for row in csv.reader(_in, delimiter="\t"):
				mge = row[3].split(";")[0]
				mgelen, matches = map(int, row[-3:-1])
				d.setdefault(mge, []).append((matches / mgelen, matches, mgelen, i))
	
	header = ["mge", "length", "longest_match", "longest_match_cov", "n_25", "n_50", "n75", "n_longest"]

	with open(sys.argv[2], "wt") as _out:
		print(*header, sep="\t", file=_out)
		for mge, alignments in sorted(d.items()):
			alignments = sorted(alignments)
			longest_cov, longest_length = alignments[-1][:2]
			n_25 = sum(1 for cov, *_, n in alignments if cov > 0.25)
			n_50 = sum(1 for cov, *_, n in alignments if cov > 0.50)
			n_75 = sum(1 for cov, *_, n in alignments if cov > 0.75)
			n_longest = sum(1 for cov, *_, n in alignments if cov == longest_cov)

			print(mge, alignments[0][-2], longest_length, longest_cov, n_25, n_50, n_75, n_longest, sep="\t", file=_out)





if __name__ == "__main__":
	main()