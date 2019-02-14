"""
This script finds the motifs in a sequence, report the locations in a bed file
and density in a bedgraph file.

Example for finding m6A motifs in the mmXist reference genome. The motif is
present once in every 17918/333 = 53.8nt. Therefore, I chose 100nt windows
python ~/Documents/scripts/motif.py DRACH \
~/Documents/chang/xist/structure/mmXist_igv/mmXist.fa \
~/Documents/chang/xist/structure/mmXist_igv/mmXist_DRACH.bed 300 50 \
~/Documents/chang/xist/structure/mmXist_igv/mmXist_DRACH_density.bedgraph

For RBP motifs
python ~/Documents/scripts/motif.py CCCC \
~/Documents/chang/xist/structure/hsXIST_igv/mmXist.fa \
~/Documents/chang/xist/structure/mmXist_igv/mmXist_HNRNPK_CCCC.bed 300 50 \
~/Documents/chang/xist/structure/mmXist_igv/mmXist_HNRNPK_CCCC_density.bedgraph

python ~/Documents/scripts/motif.py CCCC \
~/Documents/chang/xist/structure/hsXIST_igv/hsXIST.fa \
~/Documents/chang/xist/structure/hsXIST_igv/hsXIST_HNRNPK_CCCC.bed 300 50 \
~/Documents/chang/xist/structure/hsXIST_igv/hsXIST_HNRNPK_CCCC_density.bedgraph

python ~/Documents/scripts/motif.py CGG \
~/Documents/chang/xist/structure/hsXIST_igv/hsXIST.fa \
~/Documents/chang/xist/structure/hsXIST_igv/hsXIST_RBM22_CGG.bed 300 50 \
~/Documents/chang/xist/structure/hsXIST_igv/hsXIST_RBM22_CGG_density.bedgraph

python ~/Documents/scripts/motif.py GTRTG \
~/Documents/chang/xist/structure/hsXIST_igv/hsXIST.fa \
~/Documents/chang/xist/structure/hsXIST_igv/hsXIST_TARDBP_GTRTG.bed 300 50 \
~/Documents/chang/xist/structure/hsXIST_igv/hsXIST_TARDBP_GTRTG_density.bedgraph

python ~/Documents/scripts/motif.py TGT \
~/Documents/chang/xist/structure/hsXIST_igv/hsXIST.fa \
~/Documents/chang/xist/structure/hsXIST_igv/hsXIST_KHSRP_TGT.bed 300 50 \
~/Documents/chang/xist/structure/hsXIST_igv/hsXIST_KHSRP_TGT_density.bedgraph

"""

import sys, re
if len(sys.argv) < 7:
    print "Usage: python motif.py motif inputfastq outputbed wsize step density"
    sys.exit()


iupac = {"A": "A", \
         "C": "C", \
         "G": "G", \
         "T": "T", \
         "U": "U", \
         "W": "AT", \
         "S": "CG", \
         "M": "AC", \
         "K": "GT", \
         "R": "AG", \
         "Y": "CT", \
         "B": "CGT", \
         "D": "AGT", \
         "H": "ACT", \
         "V": "ACG", \
         "N": "ACGT" }
         
motif = sys.argv[1]
infile = open(sys.argv[2], 'r')
outbed = open(sys.argv[3], 'w')
wsize = int(sys.argv[4])
step = int(sys.argv[5])
outbedgraph = open(sys.argv[6], 'w')
indata = infile.readlines() #assume this is a single sequence
chrom = [line.strip('\n').strip('>') for line in indata if line[0] == '>'][0]
seq = ''.join([line.strip('\n') for line in indata if line[0] != ">"])
infile.close()
print "\nLength of input sequence:", len(seq)


#make a regex string for the motif
expandedmotif = ''
for i in motif: expandedmotif += ("[" + iupac[i] + "]")
print "Motif to search for:", expandedmotif


#output locations of the motifs in a bed file
matches = [[chrom, str(m.start(0)), str(m.end(0))] for m in \
           re.finditer(expandedmotif, seq)]
bedstring = '\n'.join(['\t'.join(i) for i in matches])
outbed.write(bedstring) 
outbed.close()



#output density of the motifs in a bedgraph file, count the center of each motif
centers = [(int(i[1]) + int(i[2])) / 2 for i in matches]
windowdict = {} #create a dictionary as we go through the centers list
bedgraphstring = ''
for x in centers:
    for i in range(wsize/step):
        start = ((x-wsize)/step+i)*step
        end = ((x-wsize)/step+i)*step + wsize
        if (start, end) in windowdict: windowdict[(start, end)] += 1
        else: windowdict[(start, end)] = 0
for key in sorted(windowdict.keys()):
    left, right = (sum(key)-step)/2, (sum(key)+step)/2
    if left >= 0 and right <= len(seq):
        record = [chrom, str(left), str(right), str(windowdict[key])]
        bedgraphstring += ('\t'.join(record) + '\n')

print
outbedgraph.write(bedgraphstring)
outbedgraph.close()
