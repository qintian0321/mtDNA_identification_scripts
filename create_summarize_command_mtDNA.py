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
exp1_sample=open("C:\Research_data\mtDNA_identification\exp1_list.txt").readlines()
exp2_sample=open("C:\Research_data\mtDNA_identification\exp2_list.txt").readlines()
experiment_group="exp2_group_result/"
data_directory="/home/cp832/qin/mtDNA_identification/"#"C:/Research_data/mtDNA_identification/"
correct_mtDNA_pos=open("C:/Research_data/mtDNA_identification/correct_mtDNAPrimer_position.csv").readlines()
output_command_list=open("C:/Research_data/mtDNA_identification/summarize_sam_commands_exp2.txt",'w')
for each_sample in exp2_sample:
    print each_sample
    each_sample=string.strip(each_sample)
    for each_group in range(0,40):
        input_samfile_name=data_directory+experiment_group+str(each_sample)+"/aligned/group"+str(each_group+1)+"_trimmed_aligned.sam"
        correct_chr_info=string.split(string.strip(correct_mtDNA_pos[each_group]),",")
        correct_chr=correct_chr_info[1]
        correct_for_pos=int(correct_chr_info[2])
        correct_rev_pos=int(correct_chr_info[3])
        output_directory=data_directory+experiment_group+(each_sample)+"/summary/"
        print >>output_command_list, "python summarize_sam_files.py"+" "+input_samfile_name+" "+str(each_group)+" "+str(correct_chr)+" "+str(correct_for_pos)+" "+str(correct_rev_pos)+" "+output_directory
output_command_list.close()
