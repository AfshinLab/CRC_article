"""
Add NAIBR orientation to BEDPE processed from LinkedSV

Usage:

	cat P19.candidates.no_type.bedpe | python add_type.py > P19.candidates.bedpe

Input should be in 7 columns format as below:

    chr1    1666975 1666976 chr1    1667159 1667160 DEL

Output would in this case be:

    Input should be in 7 columns format as below:

    chr1    1666975 1666976 chr1    1667159 1667160 DEL +-
"""

import sys


# Translate SV type to all possible breakpoint directions.
table = {"DEL": ["+-"], "INV": ["++", "--"], "DUP": ["-+"], "TRA": ["++", "--", "+-", "-+"]}

for line in sys.stdin:
	els = line.strip().split("\t")
	sv_type = els[-1]

	for naibr_type in table[sv_type]:
		print(line.strip(), naibr_type, sep="\t")
