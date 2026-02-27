#!/usr/bin/env python3

import csv
import sys

def main():
	with open(sys.argv[1], 'rt') as _in:
		for row in csv.reader(_in, delimiter='\t'):
			contig, mge_start, mge_end, mge, _, _, _, rec_start, rec_end, _, _, _, recombinase, overlap = row

			# MGE_SAMEA3888906.psa_megahit.psb_metabat2.00002_k119_43134:104944-125189;20247:18:0:0:0:19671:60;576;575M19671S
			mge, scores, alen, cigar = mge.split(";")
			_, mis, dels, ins, clips5, clips3, mapq = map(int, scores.split(":"))
			mge_start, mge_end, rec_start, rec_end, overlap, alen = map(int, (mge_start, mge_end, rec_start, rec_end, overlap, alen))
			alen = mge_end - (mge_start + 1) # mge_start comes from bed, so is 0-based!
			rlen = rec_end - rec_start + 1  # rec_start comes from gff, so is 1-based

			if alen < rlen:
				# aligned part is shorter than mge-recombinase
				continue

			if not (mge_start <= rec_start and rec_end <= mge_end):
				# mge/recombinase overlap is guaranteed by input
				# discard partial recombinase hits
				continue
			
			print(*row, sep='\t')



			

	...

if __name__ == "__main__":
	main()