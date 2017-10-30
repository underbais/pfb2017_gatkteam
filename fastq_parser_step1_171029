#!/usr/bin/env python3

import sys, os
import re

line_num_tot = ''
seq_trim = ''
qual_trim = ''

fastq_file1 = open (sys.argv[1], 'r')
fastq_file2 = open (sys.argv[2], 'r')
fastq_output1 = open('fastq_output1.txt','w')
fastq_output2 = open('fastq_output2.txt','w')


def fastq_parser(fastq_file1,fastq_file2, threshold=30):

    for line_num, line1 in enumerate(fastq_file1):      #iterating over single lines, two files
        line2 = fastq_file2.__next__()
        line_num_tot = line_num % 4
        shouldibreak = False     #setting break statement variable
        if line_num_tot == 0:
            header_line1 = line1.rstrip()   #iterating over two files
            header_line2 = line2.rstrip()
        if line_num_tot == 1:
            seq_line1 = line1.rstrip()
            seq_line2 = line2.rstrip()
        if line_num_tot == 2:
            plus_line1 = line1.rstrip()
            plus_line2 = line2.rstrip()
        if line_num_tot == 3:
            qual_line1 = line1.rstrip()
            qual_line2 = line2.rstrip()
            header_line1, seq_line_trim1, plus_line1, qual_line_trim1 = trim_function(header_line1, seq_line1, plus_line1, qual_line1)
            header_line2, seq_line_trim2, plus_line2, qual_line_trim2 = trim_function(header_line2, seq_line2, plus_line2, qual_line2)
            #executing defined trim_function
            if len(seq_line_trim1) and len(seq_line_trim2) >=threshold:
                fastq_output1.write(header_line1 + '\n' +  seq_line_trim1 + '\n' + plus_line1 + '\n' + qual_line_trim1 +'\n')
                fastq_output2.write(header_line2 + '\n' +  seq_line_trim2 + '\n' + plus_line2 + '\n' + qual_line_trim2 +'\n')
                #print('Accept')
            #if len(seq_line_trim1) and len(seq_line_trim2) <= threshold:
                #print ('not okay')

def trim_function(header_line, seq_line, plus_line, qual_line):
# Trims sequence and quality data, outputs new fastq file  

        for pos, qual_c in enumerate(qual_line):   #sliding window function
            windows = qual_line[pos:pos+5]
            score_count = 0
            shouldibreak = False     #break statement to end trim loop
            for window in windows:
                if ord(window)-33 < 20:    #criteria for score
                    score_count += 1
                    if score_count >= 3:
                        seq_trim = seq_line[0:pos]   #defining new function of trimmed lines
                        qual_trim = qual_line[0:pos]

                        return (header_line, seq_trim, plus_line, qual_trim)   #pulling information for trimmed seqs
                        shouldibreak = True    #breaking out of score loop
                        break

            if shouldibreak is True:    #breaking out of score loop
                shouldibreak = True
                break
        if shouldibreak is False:    #addressing seqs that are not trimmed
            seq_trim = seq_line
            qual_trim = qual_line

            return (header_line, seq_trim, plus_line, qual_trim)    #pulling information for untrimmed seqs

fastq_parser(fastq_file1,fastq_file2)    #
