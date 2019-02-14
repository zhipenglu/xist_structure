"""
To compare the details of the PARIS and fRIP-seq, I need to extend the fRIP-seq
mapped reads inward by the average size of each fragment (~150nt), minus the
sequenced tag (31nt on one side). The mode was determined by examining the
distribution, first for EZH2rep1 using the shortdist.pdf figure (150-31=119).
No curve fitting was done for this purpose. The extension was done using the
following commands from bedtools. A custom script frip_extend_longpairs.py was
written for this purpose.

I should use the average size in the range 0-400, not the mode in this step. 

Usage:
python ~/Documents/scripts/frip/frip_extend_longpairs.py 119 \
EZH2rep1gap_hsXIST_geometric_longanchors_LG2.sam \
EZH2rep1gap_hsXIST_geometric_longanchors_LG2extend.sam

"""

import sys, re

if len(sys.argv) < 4:
    print "Usage: python frip_extend_longpairs.py length inputsam outputsam"
    sys.exit()

length = int(sys.argv[1])
inputsam = open(sys.argv[2], 'r')
outputsam = open(sys.argv[3], 'w')
outstring = ''

for line in inputsam:
    if line[0] == '@':
        outstring += line
        continue
    sam = line.split()
    cigar = sam[5]
    cigars = re.findall('[0-9]+[MINDSHP=]', cigar)
    matches = [int(match[0:-1]) for match in cigars]
    if len(matches) ==3: # only consider simple cigar strings: M N M
        matches[0] += length
        matches[2] += length
        matches[1] -= 2*length
        newcigars = []
        for segment in cigars:
            if segment[-1] == 'M':
                newcigars.append(str(int(segment.strip('M'))+length) + 'M')
            elif segment[-1] == 'N':
                newcigars.append(str(int(segment.strip('N'))-2*length) + 'N')
            else: newcigars.append(segment)
        sam[5] = ''.join(newcigars)
        sam[9] = sam[9][:31] + 'N'*2*length + sam[9][31:]
        sam[10] = sam[10][:31] + '~'*2*length + sam[10][31:]
        outstring += '\t'.join(sam) + '\n'


outputsam.write(outstring)
inputsam.close()
outputsam.close()
