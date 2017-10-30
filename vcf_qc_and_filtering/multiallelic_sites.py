#!/usr/bin/env python3

import sys

input_vcf = sys.argv[1]
fh = open(input_vcf, "r")

AD_list = []
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
            AD_list.append(AD)

            
#print(AD_list)
multiallelic_list=[]
for element in AD_list:
    n = (element.split(','))
    if (len(n)) != 2:
        multiallelic_list.append(element)

print(len(multiallelic_list))
