# sgRNA-finder_general.py
import sys
import re

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def find_choose_sgRNA_general(seq1, minEditPos, maxEditPos, maxEditDistance, wtdelcounter, seq2, rgn):
    chosenOrientation = None
    chosenSG = None
    chosenCutPos = None
    chosenDistance = None
    gcPctg = None
    sghash = {}

    chars = list(seq1)

    # Check sense orientation
    for i in range(21, len(chars) - 1):
        tempsg = ""
        if rgn == "Cas9-NGG":
            if chars[i] == "G" and chars[i + 1] == "G":
                tempsg = "".join(chars[i - 21:i - 1])
        elif rgn == "Cas9-NG":
            if chars[i] == "G":
                tempsg = "".join(chars[i - 21:i - 1])
        elif rgn == "Cas9-SpRY":
            tempsg = "".join(chars[i - 21:i - 1])

        if tempsg and len(tempsg) == 20:
            gcCounter = sum(1 for base in tempsg if base in "GC")
            tempgcPctg = gcCounter / len(tempsg) * 100

            if rgn == "Cas9-NGG":
                match = re.search(f"{tempsg}.GG", seq1)
                tempCutPos = match.start() + 18 if match else None
            elif rgn == "Cas9-NG":
                match = re.search(f"{tempsg}.G.", seq1)
                tempCutPos = match.start() + 18 if match else None
            elif rgn == "Cas9-SpRY":
                match = re.search(f"{tempsg}...", seq1)
                tempCutPos = match.start() + 18 if match else None

            if tempCutPos is not None:
                tempDistance = minEditPos - tempCutPos
                if minEditPos >= tempCutPos and maxEditPos <= (tempCutPos + maxEditDistance) and "TTTTT" not in tempsg:
                    sghash[tempsg] = [tempsg, "sense", tempgcPctg, tempCutPos, tempDistance]

    # Check antisense orientation
    for i in range(len(chars) - 22):
        tempsg = ""
        if rgn == "Cas9-NGG":
            if chars[i] == "C" and chars[i + 1] == "C":
                tempsg = "".join(chars[i + 3:i + 23])
        elif rgn == "Cas9-NG":
            if chars[i + 1] == "C":
                tempsg = "".join(chars[i + 3:i + 23])
        elif rgn == "Cas9-SpRY":
            tempsg = "".join(chars[i + 3:i + 23])

        if tempsg and len(tempsg) == 20:
            tempsg = reverse_complement(tempsg)
            gcCounter = sum(1 for base in tempsg if base in "GC")
            tempgcPctg = gcCounter / len(tempsg) * 100

            rc_tempsg = reverse_complement(tempsg)
            if rgn == "Cas9-NGG":
                match = re.search(f"CC.{rc_tempsg}", seq1)
                tempCutPos = match.start() + 6 if match else None
            elif rgn == "Cas9-NG":
                match = re.search(f".C.{rc_tempsg}", seq1)
                tempCutPos = match.start() + 6 if match else None
            elif rgn == "Cas9-SpRY":
                match = re.search(f"...{rc_tempsg}", seq1)
                tempCutPos = match.start() + 6 if match else None

            if tempCutPos is not None:
                tempDistance = (tempCutPos - maxEditPos) + wtdelcounter
                if minEditPos >= (tempCutPos - maxEditDistance) and maxEditPos <= (tempCutPos + wtdelcounter) and "TTTTT" not in tempsg:
                    sghash[tempsg] = [tempsg, "antisense", tempgcPctg, tempCutPos, tempDistance]

    if sghash:
        for k in sghash.keys():
            orientation = sghash[k][1]
            sgRNA = k
            sgPAM = None
            rc_sgPAM = None

            if rgn == "Cas9-NGG":
                sgPAM = f"{sgRNA}.GG"
                rc_sgPAM = f"CC.{reverse_complement(sgRNA)}"
            elif rgn == "Cas9-NG":
                sgPAM = f"{sgRNA}.G."
                rc_sgPAM = f".C.{reverse_complement(sgRNA)}"
            elif rgn == "Cas9-SpRY":
                sgPAM = f"{sgRNA}..."
                rc_sgPAM = f"...{reverse_complement(sgRNA)}"

            if orientation == "sense":
                sghash[k].append(1 if re.search(sgPAM, seq2) is None else 0)
            elif orientation == "antisense":
                sghash[k].append(1 if re.search(rc_sgPAM, seq2) is None else 0)

        ranksghash = {}
        for counter, k in enumerate(sorted(sghash.keys(), key=lambda x: (sghash[x][5], sghash[x][4]), reverse=True)):
            if counter == 0:
                chosenSG = k
                chosenCutPos = sghash[k][3]
                chosenDistance = sghash[k][4]
                chosenOrientation = sghash[k][1]
                gcPctg = sghash[k][2]
            ranksghash[counter] = [k, "NA", sghash[k][3], sghash[k][1], sghash[k][4], sghash[k][5], sghash[k][2]]
            print(ranksghash[counter][3], file=sys.stderr)

        return sghash, chosenSG, chosenCutPos, chosenOrientation, chosenDistance, gcPctg, ranksghash
    else:
        return "none", "none", "none", "none", "none", "none", "none"


def find_choose_nick_sgRNA_general(seq1, minEditPos, maxEditPos, maxEditDistance, chosenCutPos, chosenOrientation, minNickDist, maxNickDist, rgn):
    chosenNickOrientation = None
    chosenNickSG = None
    chosenNickSGPos = None
    chosenNickSGDist = None
    nicksghash = {}

    chars = list(seq1)

    if chosenOrientation == "antisense":
        for i in range(21, len(chars) - 1):
            tempsg = ""
            if rgn == "Cas9-NGG":
                if chars[i] == "G" and chars[i + 1] == "G":
                    tempsg = "".join(chars[i - 21:i - 1])
            elif rgn == "Cas9-NG":
                if chars[i] == "G":
                    tempsg = "".join(chars[i - 21:i - 1])
            elif rgn == "Cas9-SpRY":
                tempsg = "".join(chars[i - 21:i - 1])

            if tempsg and len(tempsg) == 20:
                if rgn == "Cas9-NGG":
                    match = re.search(f"{tempsg}.GG", seq1)
                    tempCutPos = match.start() + 18 if match else None
                elif rgn == "Cas9-NG":
                    match = re.search(f"{tempsg}.G.", seq1)
                    tempCutPos = match.start() + 18 if match else None
                elif rgn == "Cas9-SpRY":
                    match = re.search(f"{tempsg}...", seq1)
                    tempCutPos = match.start() + 18 if match else None

                if tempCutPos is not None:
                    tempDistance = chosenCutPos - tempCutPos - 1
                    if abs(tempDistance) >= int(minNickDist) and abs(tempDistance) <= int(maxNickDist) and "TTTTT" not in tempsg:
                        nicksghash[tempsg] = [tempsg, tempCutPos, "sense", tempDistance]

    if chosenOrientation == "sense":
        for i in range(len(chars) - 22):
            tempsg = ""
            if rgn == "Cas9-NGG":
                if chars[i] == "C" and chars[i + 1] == "C":
                    tempsg = "".join(chars[i + 3:i + 23])
            elif rgn == "Cas9-NG":
                if chars[i + 1] == "C":
                    tempsg = "".join(chars[i + 3:i + 23])
            elif rgn == "Cas9-SpRY":
                tempsg = "".join(chars[i + 3:i + 23])

            if tempsg and len(tempsg) == 20:
                tempsg = reverse_complement(tempsg)
                rc_tempsg = reverse_complement(tempsg)
                if rgn == "Cas9-NGG":
                    match = re.search(f"CC.{rc_tempsg}", seq1)
                    tempCutPos = match.start() + 6 if match else None
                elif rgn == "Cas9-NG":
                    match = re.search(f".C.{rc_tempsg}", seq1)
                    tempCutPos = match.start() + 6 if match else None
                elif rgn == "Cas9-SpRY":
                    match = re.search(f"...{rc_tempsg}", seq1)
                    tempCutPos = match.start() + 6 if match else None

                if tempCutPos is not None:
                    tempDistance = tempCutPos - chosenCutPos + 1
                    if abs(tempDistance) >= minNickDist and abs(tempDistance) <= maxNickDist and "TTTTT" not in tempsg:
                        nicksghash[tempsg] = [tempsg, tempCutPos, "antisense", tempDistance]

    if nicksghash:
        for counter, k in enumerate(sorted(nicksghash.keys(), key=lambda x: abs(nicksghash[x][3] - 50))):
            if counter == 0:
                chosenNickSG = k
                chosenNickSGPos = nicksghash[k][1]
                chosenNickOrientation = nicksghash[k][2]
                chosenNickSGDist = nicksghash[k][3]
            counter += 1

        return chosenNickSG, chosenNickSGPos, chosenNickOrientation, chosenNickSGDist, nicksghash
    else:
        return "none found", "none", "none", "none", "none"


def find_pe3b_sgRNA_general(seq1, minEditPos, maxEditPos, maxEditDistance, chosenCutPos, chosenOrientation, seq2, wtdelcounter, mutdelcounter, rgn):
    chosenNickOrientation = None
    chosenNickSG = None
    chosenNickSGPos = None
    chosenNickSGDist = None
    nicksghash = {}

    chars = list(seq2)

    if chosenOrientation == "antisense":
        for i in range(21, len(chars) - 1):
            tempsg = ""
            seedsg = ""
            if rgn == "Cas9-NGG":
                if chars[i] == "G" and chars[i + 1] == "G":
                    tempsg = "".join(chars[i - 21:i - 1])
                    seedsg = "".join(chars[i - 11:i - 1])
            elif rgn == "Cas9-NG":
                if chars[i] == "G":
                    tempsg = "".join(chars[i - 21:i - 1])
                    seedsg = "".join(chars[i - 11:i - 1])
            elif rgn == "Cas9-SpRY":
                tempsg = "".join(chars[i - 21:i - 1])
                seedsg = "".join(chars[i - 11:i - 1])

            if tempsg and len(tempsg) == 20:
                if rgn == "Cas9-NGG":
                    seedPAM = f"{seedsg}.GG"
                    sgPAM = f"{tempsg}.GG"
                elif rgn == "Cas9-NG":
                    seedPAM = f"{seedsg}.G."
                    sgPAM = f"{tempsg}.G."
                elif rgn == "Cas9-SpRY":
                    seedPAM = f"{seedsg}..."
                    sgPAM = f"{tempsg}..."

                if re.search(seedPAM, seq1) is None and "TTTTT" not in tempsg:
                    match = re.search(sgPAM, seq2)
                    if match:
                        tempCutPos = match.start() + 18
                        tempDistance = chosenCutPos - tempCutPos - mutdelcounter + 1
                        nicksghash[tempsg] = [tempsg, tempCutPos, "sense", tempDistance]

    if chosenOrientation == "sense":
        for i in range(len(chars) - 22):
            tempsg = ""
            seedsg = ""
            if rgn == "Cas9-NGG":
                if chars[i] == "C" and chars[i + 1] == "C":
                    tempsg = "".join(chars[i + 3:i + 23])
                    seedsg = "".join(chars[i + 3:i + 13])
            elif rgn == "Cas9-NG":
                if chars[i + 1] == "C":
                    tempsg = "".join(chars[i + 3:i + 23])
                    seedsg = "".join(chars[i + 3:i + 13])
            elif rgn == "Cas9-SpRY":
                tempsg = "".join(chars[i + 3:i + 23])
                seedsg = "".join(chars[i + 3:i + 13])

            if tempsg and len(tempsg) == 20:
                tempsg = reverse_complement(tempsg)
                rc_tempsg = reverse_complement(tempsg)

                if rgn == "Cas9-NGG":
                    seedPAM = f"CC.{seedsg}"
                    sgPAM = f"CC.{rc_tempsg}"
                elif rgn == "Cas9-NG":
                    seedPAM = f".C.{seedsg}"
                    sgPAM = f".C.{rc_tempsg}"
                elif rgn == "Cas9-SpRY":
                    seedPAM = f"...{seedsg}"
                    sgPAM = f"...{rc_tempsg}"

                if re.search(seedPAM, seq1) is None and "TTTTT" not in tempsg:
                    match = re.search(sgPAM, seq2)
                    if match:
                        tempCutPos = match.start() + 6
                        tempDistance = tempCutPos - chosenCutPos + 1 - wtdelcounter
                        nicksghash[tempsg] = [tempsg, tempCutPos, "antisense", tempDistance]

    if nicksghash:
        for counter, k in enumerate(sorted(nicksghash.keys(), key=lambda x: nicksghash[x][3])):
            if counter == 0:
                chosenNickSG = k
                chosenNickSGPos = nicksghash[k][1]
                chosenNickOrientation = nicksghash[k][2]
                chosenNickSGDist = nicksghash[k][3]
            counter += 1

        return chosenNickSG, chosenNickSGPos, chosenNickOrientation, chosenNickSGDist, nicksghash
    else:
        return "none found", "none", "none", "none", "none"




