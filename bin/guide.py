#!/usr/bin/env python3

import re
import sys

CIGAR_RE = re.compile('r([0-9]+[MIDNSH=X])')


def parse_cigar(start, cigar):
	p = start
	mis, dels, ins = 0, 0, 0
	
	# 662S23=1X5=1X2=1X26=1X90=1X10=1X137=1X21=1X17=
	for match in CIGAR_RE.findall(cigar):
		n, op = int(match[:-1]), match[-1]
		if op in ("S", "H"):
			if p == start:
				start += n
			else:
				end = p
				length = p + n
			# p += n
		elif op in ("=", "X", "D", "N"):
			p += n
			mis += n * (op == "X")
			dels += n * (op in ("D", "N"))
		elif op == "I":
			ins += n

	return start, p, length, mis, dels, ins

		

		







def main():
	with open(sys.argv[1], "rt") as _in:
		for line in _in:
			if line[0] != "@":
				# MGE_GCA_006692365.1_593905.SAMN03466345.AAFZNY010000105:439-1439        0       k119_180526     1       60      662S23=1X5=1X2=1X26=1X90=1X10=1X137=1X21=1X17=
				fields = line.strip().split("\t")
				contig, start, cigar = fields[2], int(fields[3]), fields[5]
				mge_id, mapq = fields[0], int(fields[4])

				start, end, *scores = parse_cigar(start, cigar)
				scores = ":".join(map(str, (scores + [mapq,])))
				print(contig, start - 1, end, f"{mge_id};{scores}", sep="\t")



	


if __name__ == "__main__":
	main()