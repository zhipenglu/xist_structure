"""
pe2gap.py

Copyright: Zhipeng Lu, 2018, zhipengluchina@gmail.com
This script converts the concatenated paired-end reads to gapped reads. In
other words, the second part is reverse-complemented and quality is reversed.
The first tests were performed on fRIP-seq data from Henderickson et al. 2016. 

example:
python ~/zhipeng/bin/pe2gap.py 32 test.fastq testgap.fastq

example concatenated paired-end read:
@SRR1976644.1 HISEQ:795:H7D08ADXX:1:1101:1096:2244 length=62
GTCCCGCGGGTCTGTCTCTTGCTTCAACAGTGCAAGCGGAGAGGAANTCGGAGGGCGGTCGG
+SRR1976644.1 HISEQ:795:H7D08ADXX:1:1101:1096:2244 length=62
BBBFFFFFFFFFFIIIIIIIIIIIIIIIIIFBBBFFFFFFFFFFFI#0BFFFIIIIIIIFFF

example output gapped read:
@SRR1976644.1 HISEQ:795:H7D08ADXX:1:1101:1096:2244 length=62
GTCCCGCGGGTCTGTCTCTTGCTTCAACAGTCCGACCGCCCTCCGANTTCCTCTCCGCTTGC
+SRR1976644.1 HISEQ:795:H7D08ADXX:1:1101:1096:2244 length=62
BBBFFFFFFFFFFIIIIIIIIIIIIIIIIIFFFFIIIIIIIFFFB0#IFFFFFFFFFFFBBB
"""

import sys

if len(sys.argv) < 4:
    print "Usage: python pe2gap.py start2 pe.fastq gap.fastq"
    sys.exit()

start2 = int(sys.argv[1])
pefastq = open(sys.argv[2], 'r')
gapfastq = open(sys.argv[3], 'w')

def rc(seq):
    rule = {"A":"T", "a":"T", "T":"A", "t":"A", "U":"A", "u":"A", "G":"C", \
            "g":"C", "C":"G", "c":"G", "Y":"R", "y": "R", "R":"Y", "r":"Y", \
            "N":"N", "n":"N", "-":"-"}
    return "".join([rule[base] for base in reversed(seq)])

i = 0
outstring = ''
for line in pefastq:
    if i%4 == 0: outstring += line
    if i%4 == 1: outstring += line[0:start2-1] + rc(line[start2-1:-1]) + '\n'
    if i%4 == 2: outstring += '+\n'
    if i%4 == 3: outstring += line[0:start2-1] + line[start2-1:-1][::-1] + '\n'
    i += 1

gapfastq.write(outstring)
pefastq.close()
gapfastq.close()
