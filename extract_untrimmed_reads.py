#-------------------------------------------------------------------------------
# Name:        extract_unprocessed_reads
# Purpose:
#
# Author:      chris
#
# Created:     29/03/2018
# Copyright:   (c) chris 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import sys
import gzip
import string
read_fastq=gzip.open(sys.argv[1]).readlines()
untrimmed_id=open(sys.argv[2]).readlines()
for i in range(0,len(read_fastq),4):#for each read
    if read_fastq[i] in untrimmed_id:
        print string.strip(read_fastq[i])
        print string.strip(read_fastq[i+1])
        print string.strip(read_fastq[i+2])
        print string.strip(read_fastq[i+3])

