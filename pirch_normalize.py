"""
pirch_normalize.py

Input pirch bedgraph files, output normalized bedgraph in 100nt windows. 
Start test with mmXist in NPC cells PIRCh data, and then accomodate other data. 

This script normalizes all the PIRCh tracks using the following steps.
1. read all individual bedgraph files 
2. expand bedgraph to each nt position
3. normalize against the input controls
4. take the geometric average of each pair for display to remove the noise. 

The multibedgraph format: hsXIST 0 1 RBP1_SMinput RBP1_eCLIP1 RBP1_eCLIP2 ...
The header file is as follows for NPC PIRCh data:
NPC_H3K27Ac_rep1 NPC_H3K27Ac_rep2 NPC_H3K27Me3_rep1 NPC_H3K27Me3_rep2
NPC_H3K4Me3_rep1 NPC_H3K4Me3_rep2 NPC_IgG_rep1 NPC_IgG_rep2
NPC_Input_rep1 NPC_Input_rep2

example:
cd /Users/lu/Documents/chang/other/pirch/npc/npc_fastq
python ~/Documents/scripts/pirch/pirch_normalize.py mmXist 100

cd /Users/lu/Documents/chang/eCLIP/fripdata/
python ~/Documents/scripts/pirch/pirch_normalize.py hg38XIST 100

"""

import sys, math, os
import numpy as np

if len(sys.argv) < 2:
    print "Usage: python pirch_normalize.py RNA windowsize"
    print "RNA: mmXist"
    sys.exit()

    
RNA = sys.argv[1]
windowsize = int(sys.argv[2])
bgs = [f for f in os.listdir('./') if ("normBCD.bedgraph" in f and RNA in f)]
print bgs
inputbgs = [f for f in bgs if ("_wt_in_" in f)]
print "Input bedgraph files:", inputbgs
#set the bed file record here: 
hg38XIST = 'chrX 73820650 73852753 XIST 999 - 73820650 73852753 0 6 \
            7334,164,209,137,64,11372 0,8417,10415,12587,16789,20731'
mm10Xist = 'chrX 103460372 103483233 Xist 999 - 103460372 103483233 0 7 \
            7664,155,147,211,132,91,9518 0,8445,8927,9218,10177,10458,13343'
mmXist = 'mmXist 0 17918 mmXist 999 + 0 17918 0 1 17918 0'
hg19MALAT1 = 'chr11 65265233 65273939 MALAT1 999 + 65265233 65273939 0 1 8707 0'
hg38MALAT1 = 'chr11 65497738 65506512 MALAT1 999 + 65497738 65506512 0 1 8775 0'
mm10Malat1 = 'chr19 5795690 5802671 Malat1 999 - 5795690 5802671 0 1 6982 0'
hg38NEAT1 = 'chr11 65422798 65445538 NEAT1 999 + 65422798 65445538 0 1 22741 0'
mm10Neat1 = 'chr19 5824710 5845478 Neat1 999 - 5824710 5845478 0 1 20769 0'
rnalist = {'hg38XIST': hg38XIST, 'hg38MALAT1': hg38MALAT1, 'hg38NEAT1': \
           hg38NEAT1, 'mmXist': mmXist, 'mm10Xist': mm10Xist, 'mm10Malat1': \
           mm10Malat1, 'mm10Neat1': mm10Neat1, 'hg19MALAT1': hg19MALAT1}
print rnalist[RNA].strip('\n').split()




def expandbedgraph(bedgraph, intervals): #modified from expandbedgraph.py
    #This function expands the typical bedgraph data to intervals (like exons),
    #and output bedgraph at each nucleotide position. The bedgraph file can
    #cover more coordinates than the intervals. 
    #Bedgraph format: [[chr, start, end, value], ...], numbers are not strings
    #Interval format: [[chr, start, end, ...], ...], additional fields are OK

    bedgraphdict = {} #First read input bedgraph as a dictionary
    for record in bedgraph:
        for i in range(record[1], record[2]):
            bedgraphdict[(record[0], i)] = record[3]

    newbedgraph = [] #fill all the intervals at nt resolution, fill 0 if needed
    for interval in intervals:
        chrom = interval[0]
        for i in range(interval[1], interval[2]):
            if (chrom, i) in bedgraphdict:
                newbedgraph.append([chrom, i, i+1, bedgraphdict[(chrom, i)]])
            else: newbedgraph.append([chrom, i, i+1,0])
    return newbedgraph


def windowsum(numlist, windowsize): #output a new list by the windows
    newnumlist = []
    win = 0
    len_numlist = len(numlist)
    for i in range(len_numlist):
        if i%windowsize == 0:
            newnumlist.append(win)
            win = numlist[i]
        else:
            win += numlist[i]
    if len_numlist%windowsize != 0: newnumlist.append(win)
    return newnumlist[1:]

def bed12expand(bed12record): #convert bed12 record to exons list
    #exons: [[chr, start1, end1, cumexon1], [chr, start2, end2, cumexon2], ...]
    chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, \
           itemRgb, blockCount, blockSizes, blockStarts \
           = bed12record.strip('\n').split()
    chromStart, chromEnd = int(chromStart),int(chromEnd)
    blockCount = int(blockCount)
    blockSizes = [int(i) for i in blockSizes.split(',')]
    blockStarts = [int(i) for i in blockStarts.split(',')]
    rnalength = sum(blockSizes)
    transcript = [chrom, chromStart, chromEnd]
    exons = []
    cumexon = 0 #cummulative exon length
    for i in range(int(blockCount)):
        exonStarts = chromStart + blockStarts[i]
        exonEnds = chromStart + blockStarts[i] + blockSizes[i]
        exons.append([chrom, exonStarts, exonEnds, cumexon])
        cumexon += blockSizes[i]
    return (transcript, exons, strand, rnalength)


#take the average of controls in intervals
transcript, exons, strand, rnalength = bed12expand(rnalist[RNA])
inputavg = [0 for i in range(rnalength)]
for inputfile in inputbgs:
    inputf = open(inputfile, 'r')
    inputbg = inputf.readlines()
    inputf.close()
    inputbg = [line.strip('\n').split('\t') for line in inputbg]
    inputbg = [[x[0], int(x[1]), int(x[2]), float(x[3])] for x in inputbg if \
               len(x) ==4]
    inputbg = expandbedgraph(inputbg, exons)
    inputavg = [inputavg[i] + inputbg[i][3]/len(inputbgs) \
                for i in range(rnalength)]
inputavgwin = windowsum(inputavg, windowsize)
#for i in inputavgwin: print i
 

#read all bg files and normalize against the input controls.
for bg in bgs:
    bgfile = open(bg, 'r')
    sample = bgfile.readlines()
    bgfile.close()
    sample = [line.strip('\n').split('\t') for line in sample]
    sample = [[x[0], int(x[1]), int(x[2]), float(x[3])] for x in sample]
    sample = expandbedgraph(sample, exons) #expanded by nt position
    samplewin = windowsum([i[3] for i in sample], windowsize)
    samplewinnorm = [samplewin[i]/inputavgwin[i] if inputavgwin[i]!=0 \
                     else 0 for i in range(len(samplewin))]
    if strand == '-': samplewinnorm = samplewinnorm[::-1]

    #output to individual files
    samplebg = [[RNA if (RNA == "mmXist" or RNA == "hg38XIST") else "chr1", \
                 str(i*windowsize), str((i+1)*windowsize), \
                 str(samplewinnorm[i])] for i in range(len(samplewinnorm))]
    outstring = ('\n'.join(['\t'.join(i) for i in samplebg]) + '\n')
    outbg = open(bg.replace('.bedgraph', '') + \
                 "_norm%snt.bedgraph"%windowsize, 'w')
    outbg.write(outstring)
    outbg.close()
