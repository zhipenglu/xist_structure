#input file: sam files with or without headlines
#input range for the fRIP groups, no chromosome information needed
#Note only "MINDS" are considered from the CIGAR string
#Author: Zhipeng Lu, zhipengluchina@gmail.com
#2018-02-17

"""
python ~/Documents/scripts/frip_subset_paris.py \
AMT_Stress_trim_nodup_starhsXIST_l15p2_geometricNGmin.sam \
AMT_Stress_trim_nodup_starhsXIST_l15p2_geometricNGmin_fripsubset.sam \

example read:
#M01339    256     chr21   7814836 0       10M1I10M403206N39M93N85M18S

Post processing:
samtools view -bS -o \
AMT_Stress_trim_nodup_starhsXIST_l15p2_geometricNGmin_fripsubset.bam \
AMT_Stress_trim_nodup_starhsXIST_l15p2_geometricNGmin_fripsubset.sam
samtools sort \
AMT_Stress_trim_nodup_starhsXIST_l15p2_geometricNGmin_fripsubset.bam \
AMT_Stress_trim_nodup_starhsXIST_l15p2_geometricNGmin_fripsubset_sorted
samtools index \
AMT_Stress_trim_nodup_starhsXIST_l15p2_geometricNGmin_fripsubset_sorted.bam

"""

import sys, re
if len(sys.argv) < 3:
    print "Usage: frip_subset_paris.py input.sam output.sam"
    sys.exit()

inputsam = open(sys.argv[1], 'r')
outputsam = open(sys.argv[2], 'w')
anchor1 = [2000, 2400]
anchor2 = [11100, 11500]
anchor3 = [12900, 13300]
anchor4 = [13500, 13900]
anchor5 = [17600, 18000]
anchor6 = [18900, 19296]

output = ''
for line in inputsam:
    start, end = 0, 0
    record = line.split()

    if line[0] == "@": output += line
    else:
        start = int(record[3])
        cigars = re.findall('\d+[MIDNSPH=X]', record[5])
        readrange = 0
        for cigar in cigars:
            if cigar[-1] in ["M", "D", "N"]: readrange += int(cigar[:-1])
            elif cigar[-1] in ["I"]: readrange -= int(cigar[:-1])
            #do nothing if the cigar operation is "S" or the others [HP=X].
        end = start + readrange

    if  anchor1[0] <= start <= anchor1[1] and anchor2[0]<=end<=anchor2[1]:
        record[-1] = "DG:i:5"
        output += ('\t'.join(record) + '\n')
    elif anchor3[0] <= start <= anchor3[1] and anchor6[0] <= end <= anchor6[1]:
        record[-1] = "DG:i:4"
        output += ('\t'.join(record) + '\n')
    elif anchor4[0] <= start <= anchor4[1] and anchor5[0] <= end <= anchor5[1]:
        record[-1] = "DG:i:2"
        output += ('\t'.join(record) + '\n')
    elif anchor4[0] <= start <= anchor4[1] and anchor6[0] <= end <= anchor6[1]:
        record[-1] = "DG:i:3"
        output += ('\t'.join(record) + '\n')
    elif anchor5[0] <= start <= anchor5[1] and anchor6[0] <= end <= anchor6[1]:
        record[-1] = "DG:i:1"
        output += ('\t'.join(record) + '\n')
        
inputsam.close()
outputsam.write(output)            
outputsam.close()
