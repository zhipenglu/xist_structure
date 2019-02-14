"""
petwo2gap.py
modified from the pe2gap.py script.

This script converts the paired-end reads in two files to gapped reads in a
single file. The second read is reverse-complemented and quality is reversed.
The first tests were performed on fRIP-seq data from Henderickson et al. 2016. 


Example read in R1: 
@J00118:169:HCKLMBBXX:7:1101:24231:1279 1:N:0:NACCAG
NTTTTGACTCTTCGGGAATTTGATATTCGAGAACAGTCAAAACAAGAGAGAGTCCTTTTCCTCCTCCCGTTTTCCC
+
#AAFF<JJFJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJ

python ~/bin/petwo2gap.py pirch_1QM_1_S1_L007_R1_001.fastq \
pirch_1QM_1_S1_L007_R2_001.fastq pirch_11.fastq &

"""

import sys, itertools

if len(sys.argv) < 4:
    print "Usage: python petwo2gap.py pe1.fastq pe2.fastq gap.fastq"
    sys.exit()

pe1 = open(sys.argv[1], 'r')
pe2 = open(sys.argv[2], 'r')
gapfastq = open(sys.argv[3], 'w')

def rc(seq):
    rule = {"A":"T", "a":"T", "T":"A", "t":"A", "U":"A", "u":"A", "G":"C", \
            "g":"C", "C":"G", "c":"G", "Y":"R", "y": "R", "R":"Y", "r":"Y", \
            "N":"N", "n":"N", "-":"-"}
    return "".join([rule[base] for base in reversed(seq)])

linetotal = 0
i = 0
outstring = ''
for line1, line2 in itertools.izip(pe1, pe2):
    if i%4 == 0: outstring += line1
    if i%4 == 1: outstring += line1.strip('\n') + rc(line2.strip('\n')) + '\n'
    if i%4 == 2: outstring += '+\n'
    if i%4 == 3: outstring += line1.strip('\n') + line2.strip('\n')[::-1] + '\n'
    i += 1
    linetotal += 1
    if i == 1000000: #write content after every million lines
        gapfastq.write(outstring)
        i = 0
        outstring = ''
        print "Processed", linetotal, "lines ..."

gapfastq.write(outstring)   
pe1.close()
pe2.close()
gapfastq.close()
