"""
pca2tracks.py

This script converts the PCA analysis results for RIP/CLIP enrichment to a
minimal number of tracks for display on IGV. This approach provides more useful
information than the heatmap. The input file is *pca_array.pc.txt, and output
are the first few tracks that explain the most variance (e.g. *pc1.bedgraph).

Input format:
Interval     NAME         MEAN  PC1   PC2   ...
hsXIST_0_100 hsXIST_0_100 value value value ...

Example:
cd /Users/lu/Documents/chang/eCLIP/fripsum
python ~/Documents/scripts/pca2tracks.py \
frip_gap_hsXIST_geometric_100nt_pca_array.pc.txt 7 \
frip_gap_hsXIST_geometric_100nt_pca_array

For the PCA results from gene level, need to transpose the matrix###############
python ~/Documents/scripts/pca2tracks.py \
frip_gap_hsXIST_geometric_100nt_pca_gene.pc.txt 7 array \
frip_gap_hsXIST_geometric_100nt_pca_gene

"""

import sys

if len(sys.argv) < 4:
    print "Usage: python pca2tracks.py pca_file track_number dim output_prefix"
    print "dim: gene or array"
    sys.exit()

pcafile = open(sys.argv[1], 'r')
ntracks = int(sys.argv[2])
dimension = sys.argv[3]
outputprefix = sys.argv[4]

pcadata = pcafile.readlines()[1:] #input as a list, remove the header line
pcamatrix = [line.strip('\n').split() for line in pcadata]

meanbedgraph = open(outputprefix + "_mean.bedgraph", 'w') #output mean bedgraph
meanout = ''
for row in pcamatrix: meanout += ('\t'.join(row[0].split('_') + row[2:3]) +'\n')
meanbedgraph.write(meanout)
meanbedgraph.close()

for i in range(ntracks): #output major principal component tracks
    pctrack = open(outputprefix + '_pc' + str(i+1) + '.bedgraph', 'w')
    pctrackout = ''
    for row in pcamatrix:
        pctrackout += ('\t'.join(row[0].split('_') + row[3+i:4+i]) + '\n')
    pctrack.write(pctrackout)
    pctrack.close()

pcafile.close()
