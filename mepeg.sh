rep=5
step=150
bed=$1
nbase=200
# Edit default reference sequence fa
fa=""

if [ "$#" -ne 1 ]; then
	echo "Usage: $0 snp bed [ reference.fa ]" >&2
	exit 1
fi

mkdir -p _mepeg_tmp

# Input chr, pos, and base(s) to change in bed4 format
# Extend by n bases each way
echo Extracting wild type and mutant sequences from human genome hg 38 1>&2
awk -v n=$nbase '{print $1"\t"$2-n"\t"$2+n"\t"$4}' $bed > _mepeg_tmp/bed

bedtools getfasta -fi $fa -bed _mepeg_tmp/bed -bedOut > _mepeg_tmp/seqbed

paste <(cut -f1-3 $bed) <(cut -f4-5 _mepeg_tmp/seqbed) > _mepeg_tmp/seq2

# Generate reference and mutated sequences as tab separated file
awk -v n=$nbase '{print toupper($5)"\t"toupper(substr($5,1,n-1)) $4 toupper(substr($5,n+$3-$2+1))}' _mepeg_tmp/seq2 > _mepeg_tmp/seq

# echo ACTTCGGTTTCCTCACTGTGGCAGGGATTGGGGGAGCTCCTGAGCCCTGGGCTGGGGTGGGTCCGGAGGGGTGGGGGCTCCTGGCTGGGCTGGGTGTGCCTGTATGTGACTCAGGCTGGATGCCCAGGCCCTGAGCTGGCGGCTTTGTTCTCCCACTGGCCGGGTAGGGTGCGTGCTGCAGCCCTGAGTCACCCGGGTGTGGCCGAGCTGCTGGCGGGACACCAGAGCCATACACCGCCACACTCACAGAGATCCACACACTGGAGCTCAGGACACGGCGCCCTGAGAGCGGAGGGCCCCCCTGGTGTGGGTGCCGTCACTGGTGCCTTTCCCTCGGTGGTGCTGTGAGACATTCCCAAGAAGCTCCTTGGTGACTAAGGGCCCT ACTTCGGTTTCCTCACTGTGGCAGGGATTGGGGGAGCTCCTGAGCCCTGGGCTGGGGTGGGTCCGGAGGGGTGGGGGCTCCTGGCTGGGCTGGGTGTGCCTGTATGTGACTCAGGCTGGATGCCCAGGCCCTGAGCTGGCGGCTTTGTTCTCCCACTGGCCGGGTAGGGTGCGTGCTGCAGCCCTGACTCACCCGGGTGTGGCCGAGCTGCTGGCGGGACACCAGAGCCATACACCGCCACACTCACAGAGATCCACACACTGGAGCTCAGGACACGGCGCCCTGAGAGCGGAGGGCCCCCCTGGTGTGGGTGCCGTCACTGGTGCCTTTCCCTCGGTGGTGCTGTGAGACATTCCCAAGAAGCTCCTTGGTGACTAAGGGCCCT  > _mepeg_tmp/seq

# Run pypegfinder
echo Designing pegRNA and nick RNA sequences 1>&2
while read -r wt mut; do
	python main.py $wt $mut >> _mepeg_tmp/pypeg 2>> _mepeg_tmp/pegfinder.log
done < _mepeg_tmp/seq

# Run peglit batch
echo Designing linker sequences 1>&2
echo spacer,scaffold,template,PBS,motif > _mepeg_tmp/peglittmp.csv
while read -r sg rt pbs pe3 pe3b; do
	echo ${sg},GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC,${rt},${pbs},CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA>> _mepeg_tmp/peglittmp.csv
done < _mepeg_tmp/pypeg

peglit --batch _mepeg_tmp/peglittmp.csv --num-repeats $rep --num-steps $step > _mepeg_tmp/peglitlog

# Output sgRNA RT PBS PE3 PE3b LN
awk 'BEGIN{print "sgRNA\tRT\tPBS\tPE3\tPE3b\tLN"}NR==FNR{if(FNR>1){split($1, a, ","); l[FNR-1]=a[6]};next}{print $0"\t"l[FNR]}' _mepeg_tmp/peglittmp_linker_designs.csv _mepeg_tmp/pypeg > _mepeg_tmp/res.txt

# 
awk 'function rc(a) {
		o = ""
		for(i = length(a); i > 0; i--) o = o c[substr(a, i, 1)]
		return(o) }
	BEGIN {c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ; OFS="\t"
		print "sgF\tsgR\texlnkF\texlnkR\tonePiece\tPE3\tPE3b"}
	NR>1{sg = $1; sgF = sg; sgR = rc(sg)
		if(substr(sg,1,1)!="G") { sgF="g"sgF; sgR=sgR"c"}
		sgF="cacc"sgF"gttttaga"; sgR="tagctctaaaac"sgR
		rt = $2; pbs = $3; ss3e = rt pbs; lnk = $6
		exlnk= ss3e lnk; exlnkF = "gtgc"exlnk; exlnkR = "cgcg" rc(exlnk)
		pe3f = $4 "gttttagagctagaaatagcaag"
		pe3fb = $5 "gttttagagctagaaatagcaag"
		if(substr(pe3f,1,1)!="G") pe3f = "g"pe3f
		if(substr(pe3fb,1,1)!="G") pe3fb = "g"pe3fb
		if(substr(sg,1,1)!="G") { sg="g"sg }
		op = "gcatatGGTCTCtcacc" sg "gttttaga"
		op = op "GCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGgtgc"
		op = op exlnk "CGCGaGAGACCaatagg"
		print sgF,sgR,exlnkF,exlnkR,op,pe3f,pe3fb
	}' _mepeg_tmp/res.txt

echo Complete 1>&2
rm -r _mepeg_tmp
