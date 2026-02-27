#!/usr/bin/env python3

import csv
import sys

from collections import Counter


def main():

	# recombinase_anchors = {}

	key = None
	n_aln, coverage = 0, Counter()

	with open(sys.argv[1], 'rt') as _in, open(sys.argv[2], 'wt') as _out:
		for row in csv.reader(_in, delimiter='\t'):
			contig, mge_start, mge_end, mge, _, _, _, rec_start, rec_end, _, _, _, recombinase, overlap = row

			# MGE_SAMEA3888906.psa_megahit.psb_metabat2.00002_k119_43134:104944-125189;20247:18:0:0:0:19671:60;576;575M19671S
			mge, scores, alen, cigar = mge.split(";")
			_, mis, dels, ins, clips5, clips3, mapq = map(int, scores.split(":"))
			mge_start, mge_end, rec_start, rec_end, overlap, alen = map(int, (mge_start, mge_end, rec_start, rec_end, overlap, alen))
			alen = mge_end - (mge_start + 1) # mge_start comes from bed, so is 0-based!
			rlen = rec_end - rec_start + 1  # rec_start comes from gff, so is 1-based

			if "H" in cigar:
				# supplementary/secondary alignment
				continue

			if alen < rlen:
				# aligned part is shorter than mge-recombinase
				continue

			if not (mge_start <= rec_start and rec_end <= mge_end):
				# mge/recombinase overlap is guaranteed by input
				# discard partial recombinase hits
				continue
			
			new_key = (contig, mge_start + 1, mge_end, rec_start, rec_end, recombinase.split(";")[0].split("=")[1])
			if new_key != key:
				if key is not None:
					fr_coverage = Counter({k: v/n_aln for k, v in coverage.items()})
					for c_start in range(rec_start, min(coverage), -1):
						if fr_coverage[c_start] < 0.5:
							c_start += 1
							break
					for c_end in range(rec_end, max(coverage) + 1):
						if fr_coverage[c_end] < 0.5:
							c_end -= 1
							break
					
					rec_coverage = sum(coverage[c] for c in range(rec_start, rec_end + 1)) / n_aln
					hc_mge_coverage = sum(coverage[c] for c in range(c_start, c_end + 1)) / n_aln
					lc_mge_coverage = sum(coverage[c] for c in range(mge_start, mge_end + 1)) / n_aln

					print(
						*key, n_aln, round(hc_mge_coverage, 3), round(lc_mge_coverage, 3), round(rec_coverage, 3),
						file=_out,
						sep="\t"
					)
				key = new_key
				n_aln, coverage = 0, Counter()

			n_aln += 1
			coverage.update(range(mge_start + 1, mge_end + 1))
			# ra = recombinase_anchors.setdefault(, [0, Counter()])
			# ra[0] += 1
			# ra[1].update(range(mge_start + 1, mge_end + 1))

			# print(*row, sep='\t')
		if key is not None:
			fr_coverage = Counter({k: v/n_aln for k, v in coverage.items()})
			for c_start in range(rec_start, min(coverage), -1):
				if fr_coverage[c_start] < 0.5:
					c_start += 1
					break
			for c_end in range(rec_end, max(coverage) + 1):
				if fr_coverage[c_end] < 0.5:
					c_end -= 1
					break
			
			rec_coverage = sum(coverage[c] for c in range(rec_start, rec_end + 1)) / n_aln
			hc_mge_coverage = sum(coverage[c] for c in range(c_start, c_end + 1)) / n_aln
			lc_mge_coverage = sum(coverage[c] for c in range(mge_start, mge_end + 1)) / n_aln

			print(
				*key, n_aln, round(hc_mge_coverage, 3), round(lc_mge_coverage, 3), round(rec_coverage, 3),
				file=_out,
				sep="\t"
			)

		# with open(sys.argv[2], "wt") as _out:

		# 	for (contig, mge_start, mge_end, rec_start, rec_end, recombinase), (n_aln, coverage) in recombinase_anchors.items():
		# 		fr_coverage = Counter({k: v/n_aln for k, v in coverage.items()})
		# 		for c_start in range(rec_start, min(coverage), -1):
		# 			if fr_coverage[c_start] < 0.5:
		# 				c_start += 1
		# 				break
		# 		for c_end in range(rec_end, max(coverage) + 1):
		# 			if fr_coverage[c_end] < 0.5:
		# 				c_end -= 1
		# 				break
				
		# 		rec_coverage = sum(coverage[c] for c in range(rec_start, rec_end + 1)) / n_aln
		# 		hc_mge_coverage = sum(coverage[c] for c in range(c_start, c_end + 1)) / n_aln
		# 		lc_mge_coverage = sum(coverage[c] for c in range(mge_start, mge_end + 1)) / n_aln

		# 		print(
		# 			contig, c_start, c_end, mge_start, mge_end, rec_start, rec_end, recombinase, n_aln, round(hc_mge_coverage, 3), round(lc_mge_coverage, 3), round(rec_coverage, 3),
		# 			file=_out,
		# 			sep="\t"
		# 		)


if __name__ == "__main__":
	main()