# PBS-RT-design.py
import sys
import re

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def find_pbs(chosenSG, gcPctg):
    pbshash = {}
    pbstable = '<table style ="width:30%; float = left">'
    pbstable += "<tr><th>Length</th><th>PBS_Seq</th></tr>"

    print("Length\tPBS sequence", file=sys.stderr)
    for i in range(8, 18):
        pbsSeq = chosenSG[17 - i:17]
        pbsSeq = reverse_complement(pbsSeq)
        pbshash[i] = pbsSeq
        print(f"{i}\t{pbsSeq}", file=sys.stderr)
        pbstable += f"<tr><td>{i}</td><td>{pbsSeq}</td></tr>"

    chosenPBSlen = 13 - ((gcPctg - 55) / 5)
    chosenPBSlen = int(chosenPBSlen + 0.5)

    if chosenPBSlen < 8:
        chosenPBSlen = 8
    if chosenPBSlen > 17:
        chosenPBSlen = 17

    chosenPBS = pbshash[chosenPBSlen]
    pbstable += "</table>"
    print(f"\nChose {chosenPBS} of length {chosenPBSlen} as the PBS.\n\n", file=sys.stderr)
    return pbshash, chosenPBS, chosenPBSlen, pbstable


def find_RT(align2, chosenOrientation, minEditPos, maxEditPos, chosenCutPos, wtdelcounter):
    minimalEditLen = None
    seq2chars = list(align2)

    if chosenOrientation == "sense":
        minimalEditLen = maxEditPos - chosenCutPos + 1
    elif chosenOrientation == "antisense":
        minimalEditLen = chosenCutPos - minEditPos + 1 + wtdelcounter

    print(f"minimal edit length: {minimalEditLen}", file=sys.stderr)

    delCounter = seq2chars.count("-")
    minimalRTLen = minimalEditLen - delCounter

    rthash = {}
    allrthash = {}
    rtlengths = []
    rttable = '<table style ="width:50%; float = left">'
    rttable += "<tr><th>Length</th><th>TemplateSequence</th><th>Warnings</th>"
    rtcounter = 0

    if chosenOrientation == "sense":
        if minimalEditLen < 10:
            print(f"The edit distance is < 10 nt ({minimalEditLen} nt), with {delCounter} deletion base(s). Extracting all 10-16nt RT templates:", file=sys.stderr)
            print("Length\tSenseSequence\tTemplate_Seq\tWarnings", file=sys.stderr)
            for i in range(10, 17):
                align2Copy = align2.replace("-", "")
                templateSeq = align2Copy[chosenCutPos - 1:chosenCutPos - 1 + i]
                rcTemplate = reverse_complement(templateSeq)
                print(f"{i}\t{templateSeq}\t{rcTemplate}\t", file=sys.stderr)
                rttable += f"<tr><td>{i}</td><td>{rcTemplate}</td>"
                allrthash[i] = rcTemplate
                if templateSeq.endswith("G"):
                    print("First base of RT template is C; expected to have lower efficiency", file=sys.stderr)
                    rttable += "<td>First base is 'C'</td></tr>"
                else:
                    rtcounter += 1
                    rthash[i] = rcTemplate
                    rtlengths.append(i)
                    print("",file=sys.stderr)
                    rttable += "<td></td></tr>"

            if rtcounter == 0:
                for k in allrthash.keys():
                    rthash[k] = allrthash[k]
                    rtlengths.append(k)
        else:
            print(f"The edit distance is >= 10 nt ({minimalEditLen} nt), with {delCounter} deletion base(s). Returning all {minimalRTLen + 1}-{minimalRTLen + 8}nt RT templates:", file=sys.stderr)
            print("Length\tSenseSequence\tTemplateSequence\tWarnings", file=sys.stderr)
            for i in range(minimalRTLen + 1, minimalRTLen + 7):
                align2Copy = align2.replace("-", "")
                templateSeq = align2Copy[chosenCutPos - 1:chosenCutPos - 1 + i]
                rcTemplate = reverse_complement(templateSeq)
                print(f"{i}\t{templateSeq}\t{rcTemplate}\t", file=sys.stderr)
                rttable += f"<tr><td>{i}</td><td>{rcTemplate}</td>"
                allrthash[i] = rcTemplate
                if templateSeq.endswith("G"):
                    print("First base of RT template is C; expected to have lower efficiency", file=sys.stderr)
                    rttable += "<td>First base is 'C'</td></tr>"
                else:
                    rtcounter += 1
                    rthash[i] = rcTemplate
                    rtlengths.append(i)
                    print("", file=sys.stderr)
                    rttable += "<td></td></tr>"

            if rtcounter == 0:
                for k in allrthash.keys():
                    rthash[k] = allrthash[k]
                    rtlengths.append(k)

    elif chosenOrientation == "antisense":
        if minimalEditLen < 10:
            print(f"The edit distance is < 10 nt ({minimalEditLen} nt), with {delCounter} deletion base(s). Extracting all 10-16nt RT templates:", file=sys.stderr)
            print("Length\tSenseSequence\tTemplateSequence\tWarnings", file=sys.stderr)
            for i in range(10, 17):
                align2Copy = align2.replace("-", "")
                templateSeq = align2Copy[(chosenCutPos - delCounter + wtdelcounter) - i:(chosenCutPos - delCounter + wtdelcounter)]
                rcTemplate = templateSeq
                print(f"{i}\t{templateSeq}\t{rcTemplate}\t", file=sys.stderr)
                rttable += f"<tr><td>{i}</td><td>{rcTemplate}</td>"
                allrthash[i] = rcTemplate
                if templateSeq.startswith("C"):
                    print("First base of RT template is C; expected to have lower efficiency", file=sys.stderr)
                    rttable += "<td>First base is 'C'</td></tr>"
                else:
                    rtcounter += 1
                    rthash[i] = rcTemplate
                    rtlengths.append(i)
                    print("", file=sys.stderr)
                    rttable += "<td></td></tr>"

            if rtcounter == 0:
                for k in allrthash.keys():
                    rthash[k] = allrthash[k]
                    rtlengths.append(k)
        else:
            print(f"The edit distance is >= 10 nt ({minimalEditLen} nt), with {delCounter} deletion base(s). Extracting all {minimalRTLen + 1}-{minimalRTLen + 7}nt RT templates:", file=sys.stderr)
            print("Length\tSenseSequence\tTemplateSequence\tWarnings", file=sys.stderr)
            for i in range(minimalRTLen + 1, minimalRTLen + 7):
                align2Copy = align2.replace("-", "")
                templateSeq = align2Copy[(chosenCutPos - delCounter + wtdelcounter) - i:(chosenCutPos - delCounter + wtdelcounter)]
                rcTemplate = templateSeq
                print(f"{i}\t{templateSeq}\t{rcTemplate}\t", file=sys.stderr)
                rttable += f"<tr><td>{i}</td><td>{rcTemplate}</td>"
                allrthash[i] = rcTemplate
                if templateSeq.startswith("C"):
                    print("First base of RT template is C; expected to have lower efficiency", file=sys.stderr)
                    rttable += "<td>First base is 'C'</td></tr>"
                else:
                    rtcounter += 1
                    rthash[i] = rcTemplate
                    rtlengths.append(i)
                    print("", file=sys.stderr)
                    rttable += "<td></td></tr>"

            if rtcounter == 0:
                for k in allrthash.keys():
                    rthash[k] = allrthash[k]
                    rtlengths.append(k)

    rttable += "</table>"

    print("\nChoosing one RT template sequence...", file=sys.stderr)
    chosenRT = None

    if len(rtlengths) % 2 == 0:
        midRTindex = rtlengths[len(rtlengths) // 2]
    else:
        midRTindex = rtlengths[len(rtlengths) // 2 - 1]

    chosenRT = rthash[midRTindex]
    print(f"Of the candidate RT templates, chose {chosenRT} of length {midRTindex} (median size of candidates).", file=sys.stderr)
    return rthash, chosenRT, midRTindex, rttable

