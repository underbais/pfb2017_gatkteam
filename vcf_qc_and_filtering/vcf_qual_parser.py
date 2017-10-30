#!/usr/bin/env python3

# Create a script to evaluate quality of vcf file produced by our fastq -> vcf pipeline
# Generate the numbers required to create histograms for read depth and allelic balance from a single VCF

import sys

input_vcf = sys.argv[1]
fh = open(input_vcf, "r")

AFs_dict = {}
DPs_dict = {}
AFs_list = []
DPs_list = []
for line in fh:
    if not line.startswith('#'):
        field10 = line.split()[9]
        AD = field10.split(':')[1]
        DP = int(field10.split(':')[2])
        AD1 = float(AD.split(',')[0])
        AD2 = float(AD.split(',')[1])
        variant_AF = (AD1/(AD1+AD2))
        AFs_list.append(variant_AF)
        DPs_list.append(DP)
        
#        print("{}\t{}\t{}\t{:.4f}".format(AD,AD1,AD2,variant_AF))
# print(AFs_list,DPs_list)

        
