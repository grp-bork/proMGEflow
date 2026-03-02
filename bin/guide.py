#!/usr/bin/env python3

import csv
import sys

from collections import Counter
from contextlib import nullcontext


def process_recombinase(coverage, n_aln, rec_start, rec_end):
	fr_coverage = Counter({k: v/n_aln for k, v in coverage.items()})
	c_start, c_end = rec_start, rec_end
	mge_start, mge_end = min(coverage), max(coverage)
	for c in range(rec_start, min(coverage) + 1, -1):
		c_start = c
		if fr_coverage[c_start] < 0.5:
			c_start += 1
			break
	for c in range(rec_end, max(coverage) + 1):
		c_end = c
		if fr_coverage[c_end] < 0.5:
			c_end -= 1
			break
	
	rec_pileup = sum(fr_coverage[c] for c in range(rec_start, rec_end + 1))
	hc_mge_pileup = sum(fr_coverage[c] for c in range(c_start, c_end + 1))
	lc_mge_pileup = sum(fr_coverage[c] for c in range(mge_start, mge_end + 1))

	rec_coverage = rec_pileup / (rec_end - rec_start + 1)
	hc_mge_coverage = hc_mge_pileup / (c_end - c_start + 1)
	lc_mge_coverage = lc_mge_pileup / (mge_end - mge_start + 1)

	return (
		mge_start, mge_end,
		c_start, c_end,
		n_aln,
		round(rec_coverage, 3), round(lc_mge_coverage, 3), round(hc_mge_coverage, 3),
	)

def write_bed_line(key, res, stream=sys.stdout):
	# contig  rstart  rend    recombinase     mstart_lo       mend_lo mstart_hi       mend_hi n_aln   rcov    mcov_lo mcov_hi
	# k119_100061     20      505     phage_integrase 15      515     17      515     1       1.0     1.0     1.0
	data = [
		f"recombinase={key[3]},{key[1]}-{key[2]}",
		f"n_aln={res[4]}",
		f"hc_region={res[2]}-{res[3]},{res[7]}",
		f"lc_cov={res[6]}"
	]
	print(
		key[0], key[1] - 1, key[2], ";".join(data), sep="\t", file=stream,

	)



def main():

	# recombinase_anchors = {}

	key = None
	n_aln, coverage = 0, Counter()

	if sys.argv[1] == "-":
		_in, stream = sys.stdin, nullcontext()
	else:
		_in = stream = open(sys.argv[1], 'rt')

	with stream, open(f"{sys.argv[2]}.mge_candidates.tsv", 'wt') as _out, open(f"{sys.argv[2]}.mge_candidates.raw.tsv", "wt") as raw_out, open(f"{sys.argv[2]}.mge_candidates.bed", "wt") as bed_out:
		header = [
			"contig", "rstart", "rend", "recombinase",
			"mstart_lo", "mend_lo", "mstart_hi", "mend_hi",
			"n_aln", "rcov", "mcov_lo", "mcov_hi",
		]
		print(*header, sep="\t", file=_out)

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
			
			new_key = (contig, rec_start, rec_end, recombinase.split(";")[0].split("=")[1])
			if new_key != key:
				if key is not None:
					res = process_recombinase(coverage, n_aln, *key[1:3],)
					print(*key, *res, sep="\t", file=_out,)
					write_bed_line(key, res, stream=bed_out,)

					# fr_coverage = Counter({k: v/n_aln for k, v in coverage.items()})
					# c_start, c_end = key[1], key[2]
					# mge_start, mge_end = min(coverage), max(coverage)
					# for c in range(key[1], min(coverage) + 1, -1):
					# 	c_start = c
					# 	if fr_coverage[c_start] < 0.5:
					# 		c_start += 1
					# 		break
					# for c in range(key[2], max(coverage) + 1):
					# 	c_end = c
					# 	if fr_coverage[c_end] < 0.5:
					# 		c_end -= 1
					# 		break
					
					# rec_pileup = sum(fr_coverage[c] for c in range(rec_start, rec_end + 1))
					# hc_mge_pileup = sum(fr_coverage[c] for c in range(c_start, c_end + 1))
					# lc_mge_pileup = sum(fr_coverage[c] for c in range(mge_start, mge_end + 1))

					# rec_coverage = rec_pileup / (key[2] - key[1] + 1)
					# hc_mge_coverage = hc_mge_pileup / (c_end - c_start + 1)
					# lc_mge_coverage = lc_mge_pileup / (mge_end - mge_start + 1)

					# print(
					# 	*key,
					# 	mge_start, mge_end,
					# 	c_start, c_end,
					# 	n_aln,
					# 	round(rec_coverage, 3), round(lc_mge_coverage, 3), round(hc_mge_coverage, 3),
					# 	file=_out,
					# 	sep="\t",
					# )
				key = new_key
				n_aln, coverage = 0, Counter()

			n_aln += 1
			coverage.update(range(mge_start + 1, mge_end + 1))


			print(*row, sep='\t', file=raw_out,)

		if key is not None:
			res = process_recombinase(coverage, n_aln, *key[1:3],)
			print(*key, *res, sep="\t", file=_out,)
			write_bed_line(key, res, stream=bed_out,)
			# fr_coverage = Counter({k: v/n_aln for k, v in coverage.items()})
			# c_start, c_end = rec_start, rec_end
			# mge_start, mge_end = min(coverage), max(coverage)
			# for c in range(rec_start, min(coverage) - 1, -1):
			# 	c_start = c
			# 	if fr_coverage[c_start] < 0.5:
			# 		c_start += 1
			# 		break
			# for c in range(rec_end, max(coverage) + 1):
			# 	c_end = c
			# 	if fr_coverage[c_end] < 0.5:
			# 		c_end -= 1
			# 		break
			
			# rec_pileup = sum(fr_coverage[c] for c in range(rec_start, rec_end + 1))
			# hc_mge_pileup = sum(fr_coverage[c] for c in range(c_start, c_end + 1))
			# lc_mge_pileup = sum(fr_coverage[c] for c in range(mge_start, mge_end + 1))

			# rec_coverage = rec_pileup / (rec_end - rec_start + 1)
			# hc_mge_coverage = hc_mge_pileup / (c_end - c_start + 1)
			# lc_mge_coverage = lc_mge_pileup / (mge_end - mge_start + 1)

			# print(
			# 	*key,
			# 	mge_start, mge_end,
			# 	c_start, c_end,
			# 	n_aln,
			# 	round(rec_coverage, 3), round(lc_mge_coverage, 3), round(hc_mge_coverage, 3),
			# 	file=_out,
			# 	sep="\t",
			# )

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