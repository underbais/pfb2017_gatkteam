#!/usr/bin/

#sampling pair-end fastq files based on read counts specified by user

import sys,os,re

# taking filenames and number of reads as user input
R1 = str(sys.argv[1])
R2 = str(sys.argv[2])
num_reads = int(sys.argv[3])
R1_name = R1.split(".")[0]
R2_name = R2.split(".")[0]

#converting number of reads into line upper limit for the range in the loop
upper = num_reads*4+1

#creating file objects to work with
fastq1 = open(R1,"r")
fastq2 = open(R2,"r")

#creating output file objects
sample1 = open("{}/{}_sample_{}.fastq".format("samples",R1_name,num_reads),"w")
sample2 = open("{}/{}_sample_{}.fastq".format("samples",R2_name,num_reads),"w")

#looping through the the fastq objects and sampling reads
for count in range(1,upper):
	for line1 in fastq1:
		sample1.write(line1)
		line2 = fastq2.__next__()
		sample2.write(line2)
		break


