#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      chris
#
# Created:     05/04/2018
# Copyright:   (c) chris 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import string
import os
from collections import Counter
import sys
#exp1_sample=open("D:\Research_data\mtDNA_identification\exp1_list.txt").readlines()
exp2_sample=["NA18507-DTNA-CYCS-GAPDH_S9",
                "NA18507-ETFB-CYCS-GAPDH_S2",
                "NA18507-HJ-CYCS-GAPDH_S16",
                "NA18507-mt40n-all_S30",
                "NA18507-mt40n-CYCS-GAPDH_S23"]
#open("D:\Research_data\mtDNA_identification\exp2_list.txt").readlines()
experiment_group="exp2_group_result/"
data_directory="/home/cp832/qin/mtDNA_identification/"#"C:/Research_data/mtDNA_identification/"
output_command_list=open("C:/Research_data/mtDNA_identification/summarize_sam_commands_exp2_genomic_rest_v2.txt",'w')
for each_sample in exp2_sample:
    print each_sample
    each_sample=string.strip(each_sample)
    for each_group in range(40,47):
        input_samfile_name=data_directory+experiment_group+str(each_sample)+"/aligned/group"+str(each_group+1)+"_trimmed_aligned.sam"
        output_directory=data_directory+experiment_group+(each_sample)+"/summary/"
        print >>output_command_list, "python summarize_sam_files_genomic_primer.py"+" "+input_samfile_name+" "+str(each_group)+" "+output_directory
output_command_list.close()
