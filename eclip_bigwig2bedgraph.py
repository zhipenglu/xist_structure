"""
eclip_bigwig2bedgraph.py

This script extracts the bigWig data for each mature RNA and convert to a
single multibedgraph file on chr1:1-*.

Steps:
1. read the RNA bed12 record
2. extract the exons from bigwig to bedgraph format (bigWigToBedGraph)
2. remove introns, correct strand, and convert to chr1:1-*
2. expand bedgraph at each nt position (expandbedgraph function)
3. export to a multibedgraph format 

python ~/bin/eclip_bigwig2bedgraph.py XIST K562 . 100 \
hsXIST.bed eCLIP_K562_XIST.multibedgraph

python ~/Documents/scripts/eclip/eclip_bigwig2bedgraph.py ACTB K562 . 100 \
hsACTB.bed eCLIP_K562_ACTB.multibedgraph

eclip_bigwig2bedgraph.py

"""

import sys, os, subprocess
import numpy as np



if len(sys.argv) < 6:
    print "Usage: python %s RNA cell RBP window rnatrack multibedgraph" \
          %sys.argv[0]
    print "File name standard: K562/HepG2_RBP_0/1/2_pos/neg.bw"
    print "RNA: ACTB, MALAT1, NEAT1, XIST, etc."
    print "cell: K562 or HepG2"
    print "RBP: \".\" or a specific protein name"
    print "window: 1 or 100 nucleotides"
    print "rnatrack: example file name: hsACTB.bed, on chr1:1-*"
    sys.exit()
RNA = sys.argv[1]
cell = sys.argv[2]
RBP = sys.argv[3]
window = int(sys.argv[4])
rnatrack = open(sys.argv[5], 'w')
multibedgraph = open(sys.argv[6], 'w')



#predefined RNA track bed12 information
MALAT1 = 'chr11 65497738 65506512 MALAT1 999 + 65497738 65506512 0 1 8775 0'
NEAT1 = 'chr11 65422798 65445538 NEAT1 999 + 65422798 65445538 0 1 22741 0'
ACTB = "chr7 5527147 5530601 NM_001101|ACTB 999 - 5527748 5529657 0 6 \
744,182,439,240,129,78 0,856,1133,2013,2387,3376"
XIST = 'chrX 73820650 73852753 XIST 999 - 73820650 73852753 0 6 \
7334,164,209,137,64,11372 0,8417,10415,12587,16789,20731'
FTX = 'chrX 74028135 74293574 FTX 999 - 74028135 74293574 0 7 \
1416,245,69,101,155,147,223 0,3460,246206,252258,252795,253566,265216'
KCNQ1OT1 = 'chr11 2608328 2699998 KCNQ1OT1 999 - 2608328 2699998 0 1 91671 0'
HOTAIR1 = 'chr12 53962307 53968756 HOTAIR 999 - 53962307 53968756 0 6 \
1817,53,124,102,126,140 0,1923,3656,3968,4961,6309'
HOTAIR2 = 'chr12 53962307 53968914 HOTAIR 999 - 53962307 53968914 0 4 \
1817,120,102,298 0,3656,3968,6309'
HOTAIR3 = 'chr12 53962307 53974956 HOTAIR 999 - 53962307 53974956 0 6 \
1817,124,102,126,140,59 0,3656,3968,4961,6309,12590'
PVT1 = 'chr8 127890627 128101253 PVT1 999 + 127890627 128101253 0 8 \
371,169,301,130,113,137,204,274 0,48880,93276,98534,179532,205890,208742,210352'
SNHG1 = 'chr11 62851987 62855888 SNHG1 999 - 62851987 62855888 0 11 \
694,29,41,39,33,47,33,51,84,43,24 \
0,823,1087,1564,1813,2030,2531,2900,3145,3436,3877'
mm10Xist = 'chrX 103460372 103483233 Xist 999 - 103460372 103483233 0 7 \
            7664,155,147,211,132,91,9518 0,8445,8927,9218,10177,10458,13343'
rnadict = {"ACTB": ACTB, "MALAT1": MALAT1, "NEAT1": NEAT1, "XIST": XIST, \
           "KCNQ1OT1": KCNQ1OT1, 'HOTAIR1': HOTAIR1, 'HOTAIR2': HOTAIR2, \
           'HOTAIR3': HOTAIR3, 'FTX': FTX, 'PVT1': PVT1, 'SNHG1': SNHG1, \
           'mm10Xist': mm10Xist}



#process input RNA track bed12 definition
chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, \
       itemRgb, blockCount, blockSizes, blockStarts = \
       rnadict[RNA].strip('\n').split()
chromStart, chromEnd, blockCount = int(chromStart),int(chromEnd),int(blockCount)
thickStart, thickEnd = int(thickStart), int(thickEnd)
blockSizes = [int(i) for i in blockSizes.split(',')]
blockStarts = [int(i) for i in blockStarts.split(',')]
rnalength = sum(blockSizes)



#output a new track definitation of all the exons, including CDS
nchrom, nchromStart, nchromEnd, nstrand = 'chr1', '0', str(rnalength), '+'
nname, nscore, nitemRgb, nblockCount = name, score, itemRgb, blockCount
nblockSizes = blockSizes[::-1] if strand == '-' else blockSizes
nblockStarts = [0]
for i in range(blockCount-1): nblockStarts += [nblockStarts[-1]+nblockSizes[i]]
thickStartexon, thickEndexon = 0, 0
for i in range(blockCount):
    e1, e2 = chromStart+blockStarts[i], chromStart+blockStarts[i]+blockSizes[i]
    if thickStart >= e1 and thickStart <= e2:
        thickStartexon = [i, thickStart-e1, e2-thickStart]
    if thickEnd >= e1 and thickEnd <= e2:
        thickEndexon = [i, thickEnd-e1, e2-thickEnd]       
if strand == '-':
    thickStartexon, thickEndexon = [blockCount-1-thickEndexon[0], \
                                    thickEndexon[2], thickEndexon[1]], \
                                    [blockCount-1-thickStartexon[0], \
                                     thickStartexon[2], thickStartexon[1]]
nthickStart = nblockStarts[thickStartexon[0]] + thickStartexon[1]
nthickEnd = nblockStarts[thickEndexon[0]] + thickEndexon[1]
exonjunctions = '\n'.join(['\t'.join([nchrom, str(nblockStarts[i]), \
                                      str(nblockStarts[i]+1), 'J'+str(i)])\
                           for i in range(1, nblockCount)]) + '\n'
nblockSizes = ','.join([str(i) for i in nblockSizes])
nblockStarts = ','.join([str(i) for i in nblockStarts])
ones = ','.join(['1']*nblockCount)
newbed = '\t'.join([nchrom, str(nchromStart), str(nchromEnd), nname, nscore, \
                    nstrand, str(nthickStart), str(nthickEnd), nitemRgb, \
                    str(nblockCount), nblockSizes, nblockStarts]) + '\n'
#Output a track that for the mature transcript and exon junctions
rnatrack.write(newbed + exonjunctions)
rnatrack.close()



#make the first three columns of the multibedgraph file
exons = []
#format [[chr, start1, end1, cumexon1], [chr, start2, end2, cumexon2], ...]
cumexon = 0 #cummulative exon length
for i in range(int(blockCount)):
    exonStarts = chromStart + blockStarts[i]
    exonEnds = chromStart + blockStarts[i] + blockSizes[i]
    exons.append([chrom, exonStarts, exonEnds, cumexon])
    cumexon += blockSizes[i]
multiout = [["chr1", str(i), str(i+1)] for i in range(rnalength)]
windows = range(rnalength)
if window !=1:
    windows = range(rnalength/window+1) if rnalength%window else \
              range(rnalength/window)
    multiout = [["chr1", str(i*window), str((i+1)*window)] for i in windows]
    print "\nProcessing RNA: {} ...".format(RNA)




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
#Here is a test of the expandbedgraph function: 
#bedgraph1 = [["chr1", 2, 5, 0.1], ["chr1", 15, 20, 0.1], ["chr2", 2, 4, 0.1]]
#intervals1 = [["chr1", 1, 10], ["chr3", 1, 8]]
#print expandbedgraph(bedgraph1, intervals1)



def windowsum(numlist, windowsize): #output a new list by the windows
    newnumlist = []
    win = 0
    len_numlist = len(numlist)
    for i in range(len_numlist):
        if i%windowsize == 0:
            newnumlist.append(win)
            win = float(numlist[i])
        else:
            win += float(numlist[i])
    if len_numlist%windowsize != 0: newnumlist.append(win)
    return newnumlist[1:]



#Process all the input bw files. 
fstrand = 'pos' if strand == '+' else 'neg'
bigwigs = sorted([f for f in os.listdir('./') if \
                  cell in f and RBP in f and fstrand in f and '.bw' in f])
print "Number of bigWigs to process:", len(bigwigs)
allrnavalues = [] #store bedgraph values, each line contains data from one bw
for bigwig in bigwigs:
    #need to select the right strand for each RNA
    rnalist = []
    for exon in exons:
        cmd = ["bigWigToBedGraph", "-chrom={}".format(exon[0]), \
               "-start={}".format(exon[1]), "-end={}".format(exon[2]), bigwig, \
               "/dev/stdout"]
        exonout = subprocess.check_output(cmd)
        exonoutlist = [line.split('\t') for line in exonout.split('\n')]
        exonoutlist = [[x[0], int(x[1]), int(x[2]), x[3]] for x in exonoutlist \
                       if len(x) ==4]
        exonoutlist = expandbedgraph(exonoutlist, [exon]) #fill each nt position
        rnalist += exonoutlist

    rnavalues = []
    if strand == '+': rnavalues = [str(abs(float(x[3]))) for x in rnalist]
    else: rnavalues = [str(abs(float(x[3]))) for x in rnalist][::-1]
    if window != 1: rnavalues = [str(i) for i in windowsum(rnavalues, window)]
    allrnavalues.append(rnavalues)
    rnaout = [multiout[i] + [rnavalues[i]] for i in windows]
    rnaout = '\n'.join(['\t'.join(item) for item in rnaout]) + '\n'
    outname = bigwig[:-3] + '_' + RNA + '_' + str(window) + '.bedgraph'
    print "Processing", bigwig     
    outfile = open(outname, 'w')
    outfile.write(rnaout) #output individual files
    outfile.close()




multiout = np.concatenate((multiout, np.transpose(allrnavalues)), axis=1)
outstring = '\n'.join(['\t'.join(item)for item in multiout]) + '\n'
multibedgraph.write(outstring)
multibedgraph.close()
print #end the output by newline
