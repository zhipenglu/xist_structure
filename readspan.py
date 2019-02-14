"""
readspan.py: calculate the span of each read, and output the reads with big \
span currently, the cutoff is 500nt. The height of each figure can be further \
reduced. This script was originally written for the analysis of fRIP-seq data \
from Hendrickson et al. 2016

Copyright: Zhipeng Lu,2018, zhipengluchina@gmail.com

example command:

cd ~/Documents/chang/eCLIP/frip/
python ~/Documents/scripts/readspan.py DNMT1 1000 \
DNMT1rep1gap_hsXIST_geometric.sam \
DNMT1rep1gap_hsXIST_geometric_long.sam \
DNMT1rep1gap_hsXIST_geometric_shortdist.pdf \
DNMT1rep1gap_hsXIST_geometric_longdist.pdf \

for file in *gap_hsXIST_geometric.sam; do \
(python ~/Documents/scripts/readspan.py ${file%gap_hsXIST_geometric.sam} 1000 \
$file ${file%.sam}_long.sam \
${file%.sam}_shortdist.pdf ${file%.sam}_longdist.pdf &); done

"""

import sys, re, matplotlib
import matplotlib.pyplot as plt

if len(sys.argv) < 7:
    print "Usage: python readspan.py protein cutoff inputsam longreads \
shortdist longdist"
    sys.exit()
    
protein = sys.argv[1]
cutoff = int(sys.argv[2])
inputsam = open(sys.argv[3], 'r')
longread = open(sys.argv[4], 'w')
shortdist = sys.argv[5]
longdist = sys.argv[6]

######################### extact long fragment reads 
sizelist = [] #currently this output is not written to a file
longcount = 0
longstring = ''
for line in inputsam:
    if line[0] == "@":
        longstring += line
        continue
    if "XG:i:2" in line: continue
    sam = line.split()
    cigar = sam[5]
    cigars = re.findall('[0-9]+[MN]', cigar)
    matches = [int(match.strip("M").strip("N")) for match in cigars]
    sizelist.append(sum(matches)) #here matches include both 'M' and 'N'
    if sum(matches)>cutoff:
        longstring += line
        longcount += 1
sizestring = ' '.join([str(size) for size in sizelist])
inputsam.close()
longread.write(longstring)
longread.close()
print "\nSample:\t", protein
print "Total:\t", str(len(sizelist))
print "Long:\t", str(longcount)
print "Ratio\t", str(float(longcount)/len(sizelist))





#from itertools import groupby # display the distribution 
#for key, group in groupby(sorted(sizelist)):
#    print str(key) + "\t" + str(len(list(group)))

#########################plot distribution of short and long fragments
#parameters for fRIP-seq data analysis
#long_xticks = [0,2000,4000,6000,8000,10000]
#long_xlim = 10000

#parameters for PIRCh data analysis
long_xticks = [0,500,1000,1500,2000]
long_xlim = 2000

fig, ax = plt.subplots(figsize=(2.2,1.3))
plt.subplots_adjust(left=0.5, bottom=0.2)
n, bins, patches = plt.hist(sizelist, [x*10+0.5 for x in range(0,1000)], \
                            histtype='stepfilled', cumulative=False)
plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
ax.locator_params(nbins=5, axis='y')
ax.set_ylabel(protein + "\nfrequency", fontsize=15)
ax.set_xticks([0,200,400])
plt.xlim(0,400)
ylong = ax.get_ylim()[1]/100.0 #use /500 for MALAT1 and NEAT1, use /100 for XIST
plt.savefig(shortdist)

fig, ax = plt.subplots(figsize=(3,1.3))
plt.subplots_adjust(left=0.1, bottom=0.2)
n, bins, patches = plt.hist(sizelist, [x*10+0.5 for x in range(0,1000)], \
                            histtype='stepfilled', cumulative=False) 
plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
#ax.set_ylabel("frequency", fontsize=15)
ax.locator_params(nbins=5, axis='y')
ax.set_xticks(long_xticks)
plt.xlim(0,long_xlim)
plt.ylim(0,ylong)
plt.savefig(longdist)
