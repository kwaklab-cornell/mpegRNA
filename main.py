# main
import sys

from needleman_wunsch import *
from sgRNA_finder_general import *
from PBS_RT_design import *
from prepare_output_table import *
 
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def main(wt, edit, sgfile = "", c_sg = "", rgn = "Cas9-NGG", minNickDist = "40", maxNickDist = "150", maxSgCt = "3", pe3Bool=1):
    seq1 = wt.upper()
    seq2 = edit.upper()
    nw_out = needleman_wunsch(seq1, seq2)
    align1, align2, minEditPos, maxEditPos, trimmingStatus5p, trimmingStatus3p, wtdelcounter, alterhash, mutdelcounter = nw_out

    align1f = re.sub(r'(.{70})', r'\1<br>', align1)
    align2f = re.sub(r'(.{70})', r'\1<br>', align2)

    if trimmingStatus5p != 0 or trimmingStatus3p != 0:
        print(
            f"Wildtype and desired sequence are not aligned on the 5' or 3' ends.\n"
            f"Please trim off the hanging DNA bases and rerun pegFinder (see alignment below). \n\n"
            f"If you used sgRNA finder results (Broad or CRISPRScan), please also rerun the sgRNA finders using the revised wildtype sequence (if modified). \n\n"
            f"Wildtype:\n{align1f}\n\nEdited:\n{align2f}\n", file=sys.stderr
        )
        return

    if minEditPos <= 20 or (len(align1) - maxEditPos) <= 20:
        print(
            f"Desired edits are too close to 5' or 3' ends (see alignment below).\n"
            f"Please include longer arms flanking the desired edits (recommend >100 bp flanks).\n\n"
            f"If you used sgRNA finder results (Broad or CRISPRScan), please also rerun the sgRNA finders using the revised wildtype sequence (if modified). \n\n"
            f"Wildtype:\n{align1f}\n\nEdited:\n{align2f}\n", file=sys.stderr
        )
        return

    if not minNickDist.isdigit():
        print('Minimum nick distance is invalid, please enter a number in plaintext format.', file=sys.stderr)
        return

    if not maxNickDist.isdigit():
        print('Maximum nick distance is invalid, please enter a number in plaintext format.', file=sys.stderr)
        return

    if not maxSgCt.isdigit():
        print('Maximum # of sgRNAs to include in design table is invalid, please enter a number in plaintext format.', file=sys.stderr)
        return

    # maxEditDistance = 150
    maxEditDistance = 60
    if rgn == "Cas9-NG":
        maxEditDistance = 50
    elif rgn == "Cas9-SpRY":
        maxEditDistance = 20

    sgtable = '<table style ="width:60%">'
    nicksgtable = '<table style ="width:75%">'

    if not c_sg and not sgfile:
        sgdata = find_choose_sgRNA_general(seq1, minEditPos, maxEditPos, maxEditDistance, wtdelcounter, seq2, rgn)
        if sgdata[0] == "none":
            print('No candidate sgRNAs found.', file=sys.stderr)
            return

        sghash, chosenSG, chosenCutPos, chosenOrientation, chosenDistance, gcPctg, ranksghash = sgdata
        sgfoundStatus = 2
        sgtable += "<tr><th>sgRNA_Seq</th><th>CutPosition</th><th>Orientation</th><th>DistanceToEditStart</th><th>Seed/PAM_Disrupt</th><th>Rank</th><th>Chosen</th></tr>"
        for k in sorted(ranksghash.keys()):
            numrank = k + 1
            sgtable += f"<tr><td>{ranksghash[k][0]}</td><td>{ranksghash[k][2]}</td><td>{ranksghash[k][3]}</td><td>{ranksghash[k][4]}</td><td>{ranksghash[k][5]}</td><td>{numrank}</td>"
            if ranksghash[k][0] == chosenSG:
                sgtable += "<td>X</td></tr>"
            else:
                sgtable += "<td></td></tr>"
        sgtable += "</table>"

        if pe3Bool > 0:
            nickdata = find_choose_nick_sgRNA_general(seq1, minEditPos, maxEditPos, maxEditDistance, chosenCutPos, chosenOrientation, int(minNickDist), int(maxNickDist), rgn)
            chosenNickSG, chosenNickSGPos, chosenNickOrientation, chosenNickSGDist, nicksghash = nickdata
            pe3bnickdata = find_pe3b_sgRNA_general(seq1, minEditPos, maxEditPos, maxEditDistance, chosenCutPos, chosenOrientation, seq2, wtdelcounter, mutdelcounter, rgn)
            chosenNickSG3b, chosenNickSGPos3b, chosenNickOrientation3b, chosenNickSGDist3b, nicksghash3b = pe3bnickdata

            if chosenNickSG3b != "none found":
                nicksgtable += "<tr><th>Nicking sgRNA_Seq</th><th>CutPosition</th><th>Orientation</th><th>DistanceTo_pegRNA_nick</th><th>Type</th><th>Chosen</th></tr>"
                for k in sorted(nicksghash3b.keys(), key=lambda x: abs(nicksghash3b[x][3])):
                    nicksgtable += f"<tr><td>{k}</td><td>{nicksghash3b[k][1]}</td><td>{nicksghash3b[k][2]}</td><td>{nicksghash3b[k][3]}</td><td>PE3b</td>"
                    if k == chosenNickSG3b:
                        nicksgtable += "<td>X (PE3b)</td></tr>"
                    else:
                        nicksgtable += "<td></td></tr>"

            if chosenNickSG != "none found":
                if chosenNickSG3b == "none found":
                    nicksgtable += "<tr><th>Nicking sgRNA_Seq</th><th>CutPosition</th><th>Orientation</th><th>DistanceTo_pegRNA_nick</th><th>Type</th><th>Chosen</th></tr>"
                for k in sorted(nicksghash.keys(), key=lambda x: abs(nicksghash[x][3] - 50)):
                    nicksgtable += f"<tr><td>{k}</td><td>{nicksghash[k][1]}</td><td>{nicksghash[k][2]}</td><td>{nicksghash[k][3]}</td><td>PE3</td>"
                    if k == chosenNickSG:
                        nicksgtable += "<td>X (PE3)</td></tr>"
                    else:
                        nicksgtable += "<td></td></tr>"

            elif chosenNickSG == "none found" and chosenNickSG3b == "none found":
                nicksgtable += '<tr><br>No suitable secondary nicking sgRNA found</tr>'
            nicksgtable += "</table>"

    if chosenSG and sgfoundStatus == 2:
        rtdata = find_RT(align2, chosenOrientation, minEditPos, maxEditPos, chosenCutPos, wtdelcounter)
        rthash, chosenRT, chosenRTlen, rttable = rtdata

        pbsdata = find_pbs(chosenSG, gcPctg)
        pbshash, chosenPBS, chosenPBSlen, pbstable = pbsdata

        extension = chosenRT + chosenPBS
        scaffold = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"

        if c_sg:
            tabledata = prep_table_chosen(rthash, pbshash, chosenRTlen, chosenPBSlen, chosenSG, maxSgCt, rgn, chosenOrientation, gcPctg, chosenDisrupt)
        else:
            tabledata = prep_table_multi(ranksghash, align2, minEditPos, maxEditPos, wtdelcounter, maxSgCt, rgn)

        fileout = tabledata[0]

        pegRNA = chosenSG + scaffold + chosenRT + chosenPBS
        oligotable = '<table style ="width:80%; float = left">'
        oligotable += "<tr><th>OligoName</th><th>Sequence</th><th>Description</th></tr>"

        if chosenSG.startswith("G"):
            oligotable += f'<tr><td>sgF</td><td>cacc{chosenSG}gttttaga</td><td>sgRNA, forward</tr>'
            oligotable += f'<tr><td>sgR</td><td>tagctctaaaac{reverse_complement(chosenSG)}</td><td>sgRNA, reverse</td></tr>'
        else:
            oligotable += f'<tr><td>sgF</td><td>caccg{chosenSG}gttttaga</td><td>sgRNA, forward</tr>'
            oligotable += f'<tr><td>sgR</td><td>tagctctaaaac{reverse_complement(chosenSG)}c</td><td>sgRNA, reverse</td></tr>'

        oligotable += '<tr><td>scaffF</td><td>GCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCG</td><td>Scaffold, forward (invariant)</td></tr>'
        oligotable += '<tr><td>scaffR</td><td>GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTC</td><td>Scaffold, reverse (invariant)</td></tr>'
        oligotable += f'<tr><td>extensF</td><td>gtgc{extension}</td><td>3\' extension, forward</td></tr>'
        oligotable += f'<tr><td>extensR</td><td>aaaa{reverse_complement(extension)}</td><td>3\' extension, reverse</td></tr>'

        if pe3Bool > 0 and (chosenNickSG or chosenNickSG3b):
            if chosenNickSG != "none found":
                if chosenNickSG.startswith("G"):
                    oligotable += f'<tr><td>PE3_sgF</td><td>cacc{chosenNickSG}</td><td>PE3 nick sgRNA, forward</td></tr>'
                    oligotable += f'<tr><td>PE3_sgR</td><td>aaac{reverse_complement(chosenNickSG)}</td><td>PE3 nick sgRNA, reverse</td></tr>'
                else:
                    oligotable += f'<tr><td>PE3_sgF</td><td>caccg{chosenNickSG}</td><td>PE3 nick sgRNA, forward</td></tr>'
                    oligotable += f'<tr><td>PE3_sgR</td><td>aaac{reverse_complement(chosenNickSG)}c</td><td>PE3 nick sgRNA, reverse</td></tr>'

            if chosenNickSG3b != "none found":
                if chosenNickSG3b.startswith("G"):
                    oligotable += f'<tr><td>PE3b_sgF</td><td>cacc{chosenNickSG3b}</td><td>PE3b nick sgRNA, forward</td></tr>'
                    oligotable += f'<tr><td>PE3b_sgR</td><td>aaac{reverse_complement(chosenNickSG3b)}</td><td>PE3b nick sgRNA, reverse</td></tr>'
                else:
                    oligotable += f'<tr><td>PE3b_sgF</td><td>caccg{chosenNickSG3b}</td><td>PE3b nick sgRNA, forward</td></tr>'
                    oligotable += f'<tr><td>PE3b_sgR</td><td>aaac{reverse_complement(chosenNickSG3b)}c</td><td>PE3b nick sgRNA, reverse</td></tr>'

        oligotable += '</table>'
    # Run peglit to find linker sequence
    # linkers = peglit.pegLIT(chosenSG, "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC", chosenRT, chosenPBS, "CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA")

    # Print output: sgRNA, RT, PBS, PE3, PE3b
    if chosenNickSG == "none found":
        chosenNickSG = "NNNNNN"
    if chosenNickSG3b == "none found":
        chosenNickSG3b = "NNNNNN"

    print(chosenSG, chosenRT, chosenPBS, chosenNickSG, chosenNickSG3b, sep = "\t")

n = len(sys.argv)
if n > 2:
    main(sys.argv[1], sys.argv[2])


