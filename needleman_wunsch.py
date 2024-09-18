import re
import sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def needleman_wunsch(seq1, seq2):
    alignments = pairwise2.align.globalms(seq1, seq2, 1, -1, -2, -0.1)
    align1, align2, score, start, end = alignments[0]

    print(align1, file=sys.stderr)
    print(align2, file=sys.stderr)

    trimmingStatus5p = 0
    if align1[0] == "-":
        trimmingStatus5p = 1
    elif align2[0] == "-":
        trimmingStatus5p = 2

    trimmingStatus3p = 0
    if align1[-1] == "-":
        trimmingStatus3p = 1
    elif align2[-1] == "-":
        trimmingStatus3p = 2

    wtdelcounter = align1.count("-")
    mutdelcounter = align2.count("-")
    alterhash = {}

    for i in range(len(align1)):
        if align1[i] != align2[i]:
            print(f"{i+1}\t{align1[i]}\t{align2[i]}", file=sys.stderr)
            alterhash[i+1] = f"{align1[i]}>{align2[i]}"

    print(f"wildtype deletion: {wtdelcounter}", file=sys.stderr)

    minEditPos = min(alterhash.keys(), default=None)
    maxEditPos = max(alterhash.keys(), default=None)

    return align1, align2, minEditPos, maxEditPos, trimmingStatus5p, trimmingStatus3p, wtdelcounter, alterhash, mutdelcounter


