#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      chris
#
# Created:     10/04/2018
# Copyright:   (c) chris 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os
import string
exp2_file=open("C:\Research_data\mtDNA_identification\exp2_list.txt").readlines()
for each_file in exp2_file:
    each_file=string.strip(each_file)
    #ETFB_path="C:/Research_data/mtDNA_identification/exp2_group_result/"+each_file+"/summary/group41_raw_summary.txt"
    #DTNA_path="C:/Research_data/mtDNA_identification/exp2_group_result/"+each_file+"/summary/group42_raw_summary.txt"
    #HJ_path="C:/Research_data/mtDNA_identification/exp2_group_result/"+each_file+"/summary/group43_raw_summary.txt"
##    ETFB_exist=os.path.isfile(ETFB_path)
##    DTNA_exist=os.path.isfile(DTNA_path)
##    HJ_exist=os.path.isfile(HJ_path)
##
##    if ETFB_exist:
##        print each_file+" ETFB:1"
##    else:
##        print each_file+" ETFB:0"
##    if DTNA_exist:
##        print each_file+" DTNA:1"
##    else:
##        print each_file+" DTNA:0"
##    if HJ_exist:
##        print each_file+" HJ:1"
##    else:
##        print each_file+" HJ:0"
    GAPDH_path="C:/Research_data/mtDNA_identification/exp2_group_result/"+each_file+"/summary/group46_raw_summary.txt"
    CYCS_path="C:/Research_data/mtDNA_identification/exp2_group_result/"+each_file+"/summary/group47_raw_summary.txt"
    if not os.path.isfile(GAPDH_path):
        print each_file+" GAPDH: raw summary does not exist"
    if not os.path.isfile(CYCS_path):
        print each_file+" CYCS: raw summary not exist"

