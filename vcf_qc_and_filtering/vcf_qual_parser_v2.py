#!/usr/bin/env python3

# Create a script to evaluate quality of vcf file produced by our fastq -> vcf pipeline
# Generate the numbers required to create histograms for read depth and allelic balance from a single VCF
# Generate two lists (variant allele frequencies and variant read depth) for unfiltered VCF, then filter and generate the two lists again.

# Filtering parameters:
# 1) GQ >= 30
# 2) AD1 >= 3 AND AD2 >=3
# 3) DP <= 100

import sys
import statistics

input_vcf = sys.argv[1]
fh = open(input_vcf, "r")

AFs_list = []
DPs_list = []
GQ_filtered_AFs_list = []
GQ_filtered_DPs_list = []
AD_filtered_AFs_list = []
AD_filtered_DPs_list = []
DP_filtered_AFs_list = []
DP_filtered_DPs_list = []
threeX_filtered_AFs_list = []
threeX_filtered_DPs_list = []

for line in fh:
    if not line.startswith('#'):
        line = line.rstrip()
        field10 = line.split("\t")[9]
        elements_in_field10 = len(field10.split(":"))
        if elements_in_field10 == 5:
            AD = field10.split(':')[1]
            DP = int(field10.split(':')[2])
            GQ = int(field10.split(':')[3])
            AD1 = float(AD.split(',')[0])
            AD2 = float(AD.split(',')[1])
            variant_AF = (AD1/(AD1+AD2))
            AFs_list.append(variant_AF)
            DPs_list.append(DP)
            if GQ >= 30:
                GQ_filtered_AFs_list.append(variant_AF)
            if GQ >= 30:
                GQ_filtered_DPs_list.append(DP)
            if AD1 >= 3 and AD2 >= 3:
                AD_filtered_AFs_list.append(variant_AF)
            if AD1 >= 3 and AD2 >= 3:
                AD_filtered_DPs_list.append(DP)
            if DP <= 100:
                DP_filtered_AFs_list.append(variant_AF)
            if DP <= 100:
                DP_filtered_DPs_list.append(DP)
            if GQ >= 30 and AD1 >= 3 and AD2 >= 3 and DP <= 100:
                threeX_filtered_AFs_list.append(variant_AF)
            if GQ >= 30 and AD1 >= 3 and AD2 >= 3 and DP <= 100:
                threeX_filtered_DPs_list.append(DP)

                
#        else:
#            print(line)
            
#        print("{}\t{}\t{}\t{:.4f}".format(AD,AD1,AD2,variant_AF))
# print(AFs_list,DPs_list)
print(len(DPs_list))
print(len(AFs_list))  
print(len(GQ_filtered_AFs_list))
print(len(GQ_filtered_DPs_list))
print(len(AD_filtered_AFs_list))
print(len(AD_filtered_DPs_list))
print(len(DP_filtered_AFs_list))
print(len(DP_filtered_DPs_list))
print(len(threeX_filtered_AFs_list))
print(len(threeX_filtered_DPs_list))

# print(statistics.median(DPs_list))

