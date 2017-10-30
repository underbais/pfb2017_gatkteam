#!/usr/bin/env python3

# Create a script to evaluate quality of vcf file produced by our fastq -> vcf pipeline
# Generate the numbers required to create histograms for read depth and allelic balance from a single VCF
# Generate two lists (variant allele frequencies and variant read depth) for unfiltered VCF, then filter and generate the two lists again for each of the following parameters:
# Will create 10 lists total.

# Filtering parameters:
# 1) GQ >= 30
# 2) AD1 >= 3 AND AD2 >=3
# 3) DP <= 100
# 4) GQ >= 30 AND AD1 >= 3 AND AD2 >=3 AND DP <= 100

import sys
import os

def gen_QC_lists(fh):
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
    QC_list_of_lists = []
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
                if AD1+AD2 != 0:
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
    QC_list_of_lists.append(AFs_list)
    QC_list_of_lists.append(DPs_list)
    QC_list_of_lists.append(GQ_filtered_AFs_list)
    QC_list_of_lists.append(GQ_filtered_DPs_list)
    QC_list_of_lists.append(AD_filtered_AFs_list)
    QC_list_of_lists.append(AD_filtered_DPs_list)
    QC_list_of_lists.append(DP_filtered_AFs_list)
    QC_list_of_lists.append(DP_filtered_DPs_list)
    QC_list_of_lists.append(threeX_filtered_AFs_list)
    QC_list_of_lists.append(threeX_filtered_DPs_list)
    return(QC_list_of_lists)

def gen_filtered_vcf(fn, fh):
    sample_name = (os.path.basename(fn)).split(".")[0]
    filtered_vcf = open(sample_name+".filtered.vcf", "w")
    for line in fh:
        if line.startswith('#'):
            filtered_vcf.write(line)
        else:
            line = line.rstrip()
            field10 = line.split("\t")[9]
            elements_in_field10 = len(field10.split(":"))
            if elements_in_field10 == 5:
                AD = field10.split(':')[1]
                DP = int(field10.split(':')[2])
                GQ = int(field10.split(':')[3])
                AD1 = float(AD.split(',')[0])
                AD2 = float(AD.split(',')[1])
                if GQ >= 30 and AD1 >= 3 and AD2 >= 3 and DP <= 100:
                    filtered_vcf.write(line+"\n")
    return filtered_vcf
    

input_vcf = sys.argv[1]
fh = open(input_vcf, "r")

QC = gen_QC_lists(fh)
for list in QC:
    print(len(list))

fh.seek(0)    
gen_filtered_vcf(input_vcf,fh)    

fh.close()

