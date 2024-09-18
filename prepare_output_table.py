# prepare-output-table.py
import sys
from PBS_RT_design import *
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def prep_table_multi(ranksghash, align2, minEditPos, maxEditPos, wtdelcounter, maxSgCt, rgn):
    fileout = "RT_len\tRT_seq\tRT_picked\tPBS_len\tPBS_seq\tPBS_picked\t3'_extension_seq\t3'_extension_picked\textensF_oligo\textensR_oligo\tsgRNA_seq\tsgRNA_rank\tsgF_oligo\tsgR_oligo\tsg_Orientation\tsg_Seed/PAM_disrupt\tsg_GC\%\tsg_OnTargetScore\tEnzyme\n"

    for i in sorted(ranksghash.keys()): #range(int(maxSgCt)):
        chosenSG = ranksghash[i][0]
        chosenOrientation = ranksghash[i][3]
        chosenCutPos = ranksghash[i][2]
        gcPctg = ranksghash[i][6]

        rthash, chosenRT, chosenRTlen, rttable = find_RT(align2, chosenOrientation, minEditPos, maxEditPos, chosenCutPos, wtdelcounter)
        pbshash, chosenPBS, chosenPBSlen, pbstable = find_pbs(chosenSG, gcPctg)

        for k in sorted(rthash.keys()):
            for j in sorted(pbshash.keys()):
                fileout += f"{k}\t{rthash[k]}\t"
                fileout += "1\t" if k == chosenRTlen else "0\t"
                fileout += f"{j}\t{pbshash[j]}\t"
                fileout += "1\t" if j == chosenPBSlen else "0\t"

                ext = rthash[k] + pbshash[j]
                rc_ext = reverse_complement(ext)
                rank = i + 1
                fileout += f"{ext}\t"
                fileout += "1" if k == chosenRTlen and j == chosenPBSlen else "0"
                fileout += f"\tgtgc{ext}\taaaa{rc_ext}\t{chosenSG}\t{rank}\t"

                if chosenSG.startswith("G"):
                    fileout += f"cacc{chosenSG}gttttaga\t"
                    fileout += f"tagctctaaaac{reverse_complement(chosenSG)}\t"
                else:
                    fileout += f"caccg{chosenSG}gttttaga\t"
                    fileout += f"tagctctaaaac{reverse_complement(chosenSG)}c\t"

                fileout += f"{chosenOrientation}\t{ranksghash[i][5]}\t{ranksghash[i][6]}\t{ranksghash[i][1]}\t{rgn}\n"

    return fileout, 1

def prep_table_chosen(rthash, pbshash, chosenRTlen, chosenPBSlen, chosenSG, maxSgCt, rgn, chosenOrientation, gcPctg, chosenDisrupt):
    fileout = "RT_len\tRT_seq\tRT_picked\tPBS_len\tPBS_seq\tPBS_picked\t3'_extension_seq\t3'_extension_picked\textensF_oligo\textensR_oligo\tsgRNA_seq\tsgRNA_rank\tsgF_oligo\tsgR_oligo\tsg_Orientation\tsg_Seed/PAM_disrupt\tsg_GC\%\tEnzyme\n"

    for k in sorted(rthash.keys()):
        for j in sorted(pbshash.keys()):
            fileout += f"{k}\t{rthash[k]}\t"
            fileout += "1\t" if k == chosenRTlen else "0\t"
            fileout += f"{j}\t{pbshash[j]}\t"
            fileout += "1\t" if j == chosenPBSlen else "0\t"

            ext = rthash[k] + pbshash[j]
            rc_ext = reverse_complement(ext)
            rank = "preselected"
            fileout += f"{ext}\t"
            fileout += "1" if k == chosenRTlen and j == chosenPBSlen else "0"
            fileout += f"\tgtgc{ext}\taaaa{rc_ext}\t{chosenSG}\t{rank}\t"

            if chosenSG.startswith("G"):
                fileout += f"cacc{chosenSG}gttttaga\t"
                fileout += f"tagctctaaaac{reverse_complement(chosenSG)}\t"
            else:
                fileout += f"caccg{chosenSG}gttttaga\t"
                fileout += f"tagctctaaaac{reverse_complement(chosenSG)}c\t"

            fileout += f"{chosenOrientation}\t{chosenDisrupt}\t{gcPctg}\t{rgn}\n"

    return fileout, 1

