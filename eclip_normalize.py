"""
eclip_normalize.py
This script normalizes all the eCLIP tracks using the following steps.
1. Take the average of all controls on a nucleotide level
2. For each RBP, make a custom control from the average, scaled by median.
3. Normalize against the custom made control profile.
4. Cut out zero values from the output
5. The normalized output are arbitrarily placed on chr1, starting at position 1

The multibedgraph format: hsXIST 0 1 RBP1_SMinput RBP1_eCLIP1 RBP1_eCLIP2 ...
The header files are 363 and 309 elements for K562 and HepG2 data

example:
python ~/Documents/scripts/eclip/eclip_normalize.py hsXIST 100 \
header_K562_363.txt eCLIP_K562_hsXIST_all121_100nt.multibedgraph

python ~/Documents/scripts/eclip/eclip_normalize.py hg38MALAT1 100 \
header_K562_363.txt eCLIP_K562_hg38MALAT1_all121_100nt.multibedgraph

python ~/Documents/scripts/eclip/eclip_normalize.py hg38NEAT1 100 \
header_K562_363.txt eCLIP_K562_hg38NEAT1_all121_100nt.multibedgraph

python ~/Documents/scripts/eclip/eclip_normalize.py hg38MALAT1 100 \
header_HepG2_309.txt eCLIP_HepG2_hg38MALAT1_all103_100nt.multibedgraph

python ~/Documents/scripts/eclip/eclip_normalize.py hg38NEAT1 100 \
header_HepG2_309.txt eCLIP_HepG2_hg38NEAT1_all103_100nt.multibedgraph


Then use the pc2track.py command to make tracks for visualization.
python ~/Documents/scripts/pca2tracks.py \
eCLIP_HepG2_hg38MALAT1_all103_100nt_pca_array.pc.txt 7 array \
eCLIP_HepG2_hg38MALAT1_all103_100nt_pca_array

python ~/Documents/scripts/pca2tracks.py \
eCLIP_HepG2_hg38NEAT1_all103_100nt_pca_array.pc.txt 7 array \
eCLIP_HepG2_hg38NEAT1_all103_100nt_pca_array

python ~/Documents/scripts/pca2tracks.py \
eCLIP_K562_hsXIST_all121_100nt_pca_array.pc.txt 7 array \
eCLIP_K562_hsXIST_all121_100nt_pca_array

python ~/Documents/scripts/pca2tracks.py \
eCLIP_K562_hg38MALAT1_all121_100nt_pca_array.pc.txt 7 array \
eCLIP_K562_hg38MALAT1_all121_100nt_pca_array

python ~/Documents/scripts/pca2tracks.py \
eCLIP_K562_hg38NEAT1_all121_100nt_pca_array.pc.txt 7 array \
eCLIP_K562_hg38BEAT1_all121_100nt_pca_array

"""

import sys, math
import numpy as np

if len(sys.argv) < 5:
    print "Usage: python eclip_normalize.py RNA window header multibedgraph"
    print "RNA: hsXIST, hg38MALAT1 or hg38NEAT1"
    sys.exit()
RNA = sys.argv[1]
window = int(sys.argv[2])
header = open(sys.argv[3], 'r')
multibedgraph = open(sys.argv[4], 'r')
normmatrixfile = open(sys.argv[4][:-13] + 'normmatrix', 'w')



#read and process the input header and multibedgraph files
samplenames = header.readline().split() #format: K562_AARS_eCLIP1
bgdata = multibedgraph.readlines() #read in entire file as a list
bgmatrix = [line.strip('\n').split()[3:] for line in bgdata]
bgmatrix = [[0 if x=='.' else x for x in row] for row in bgmatrix]
bgmatrix = [[float(i) for i in row] for row in bgmatrix] #each row is a position
bgmatrixtrans = np.transpose(bgmatrix) #each row is a sample (Input CLIP CLIP)
nsamples = len(bgmatrix[0])
print "Number of files to process:", nsamples
npositions = len(bgmatrixtrans[0])


print "Calculating the average control for all files (K562:121, HepG2:103)"
ctrlavg = []
for row in bgmatrix: 
    avg = sum([row[i*3] for i in range(nsamples/3)])*3/nsamples
    ctrlavg.append(avg)
ctrlmedian = np.median(ctrlavg)



print "Normalizing all the eCLIP data against input ..."
eclipmatrix = [] #each row is a normalized eclip sample, no controls
RBPctrl = bgmatrixtrans[0]
lowinputlist = []
for i in range(nsamples):
    if i%3 == 0:
        RBPctrl = bgmatrixtrans[i]
        continue
    if np.median(RBPctrl) == 0:
        print "Low input (median=0)", samplenames[i]
        lowinputlist.append(samplenames[i])
        continue
    adjust = ctrlmedian/np.median(RBPctrl)
    RBPeclip = bgmatrixtrans[i] #normalize and convert "nan" and "inf" to 0
    RBPeclip = [adjust * RBPeclip[j]/ctrlavg[j] for j in range(npositions)]
    RBPeclip = [0 if (math.isinf(x)or math.isnan(x)) else x for x in RBPeclip]
    #bring down the few big numbers since these are caused by the normalization
    cutoff = sorted(RBPeclip, reverse=True)[npositions/100]
    RBPeclip = [x if x < cutoff else cutoff for x in RBPeclip]

    RBPout = '' #output individual files, see notebook 2016-09-28 
    for j in range(npositions):
        if RBPeclip[j] == 0: continue
        RBPout += ('\t'.join(["chr1",str(j*window),str((j+1)*window),\
                              str(RBPeclip[j])]) +'\n')
    filename = samplenames[i][:-6]+ RNA +'_norm_' + str(i%3) + '.bedgraph'
    outfile = open(filename, 'w')
    outfile.write(RBPout)
    outfile.close()
    eclipmatrix.append(RBPeclip) #ech row is a sample



print "Low input list", lowinputlist
#output the normalized matrix as one file
#eclipnames: 243 for K562, 207 for HepG2
eclipnames = [['Interval'] + [name for name in samplenames if \
                              ("eCLIP" in name and name not in lowinputlist)]]
positions = [['_'.join(line.split()[0:3]) for line in bgdata]]
outmatrix = eclipnames + list(np.transpose(positions + eclipmatrix))
outstring = '\n'.join(['\t'.join(row) for row in outmatrix])



print "Outputing normalized matrix file for subsequent cluster and PCA analysis"
normmatrixfile.write(outstring)
normmatrixfile.close()
   
header.close()
multibedgraph.close()





