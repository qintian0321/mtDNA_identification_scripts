#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Wenyi Qin
#
# Created:     08/04/2018
# Copyright:   (c) Wenyi Qin 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import string
import os
import sys
import timeit
start = timeit.default_timer()
input_samfile_name=sys.argv[1]
each_group=int(sys.argv[2])
pairend_flag_set=["65,129","81,161","83,163","97,145","99,147"]#all flags mean that at least the pair end read is paired and mapped to forward and reverse strain
input_samfile=open(input_samfile_name,"r").readlines()
sam_content=[]#Store every information line in the sam file
for each_line in input_samfile:
    if each_line[0]!="@":
        sam_content.append(each_line)
#Extract unique names of all paired end reads from this sample
#if there are any aligned paired end reads in this sam file
if len(sam_content)>0:
    pairend_name_all=[]
    for j in range(0,len(sam_content)):
        read_name=string.split(sam_content[j],"\t")[0]
        pairend_name_all.append(read_name)
    pairend_name_uniq=list(set(pairend_name_all))#unique paired end read names and store them into a list
    total_group_reads=len(pairend_name_uniq)#total number of paired end reads with this primer in this sample
    sam_summary_all=[]#put the summary of a paired end read in this list
    #for each paired end read of one primer group in one sample, the processing of reads starts here
    output_dict={'chr1;SameChrMapping;Unique':0,
        'chr2;SameChrMapping;Unique':0,
        'chr3;SameChrMapping;Unique':0,
        'chr4;SameChrMapping;Unique':0,
        'chr5;SameChrMapping;Unique':0,
        'chr6;SameChrMapping;Unique':0,
        'chr7;SameChrMapping;Unique':0,
        'chr8;SameChrMapping;Unique':0,
        'chr9;SameChrMapping;Unique':0,
        'chr10;SameChrMapping;Unique':0,
        'chr11;SameChrMapping;Unique':0,
        'chr12;SameChrMapping;Unique':0,
        'chr13;SameChrMapping;Unique':0,
        'chr14;SameChrMapping;Unique':0,
        'chr15;SameChrMapping;Unique':0,
        'chr16;SameChrMapping;Unique':0,
        'chr17;SameChrMapping;Unique':0,
        'chr18;SameChrMapping;Unique':0,
        'chr19;SameChrMapping;Unique':0,
        'chr20;SameChrMapping;Unique':0,
        'chr21;SameChrMapping;Unique':0,
        'chr22;SameChrMapping;Unique':0,
        'chrX;SameChrMapping;Unique':0,
        'chrY;SameChrMapping;Unique':0,
        'chrM;SameChrMapping;Unique':0
    }
    for each_name in pairend_name_uniq:
        #Time consuming step, consider improvement later
        #print each_name
        pairend_index=[k for k, x in enumerate(pairend_name_all) if x==each_name]#sometimes there would be three lines for one paired end read, suggesting supplementary alignment
        #store flags of this read into temp_flag
        temp_flag=[]
        new_index=[]
        for each_index in pairend_index:
            read_flag=string.split(sam_content[each_index],"\t")[1]
            if (int(read_flag)<2000):#do not consider supplementary alignment
                temp_flag.append(string.split(sam_content[each_index],"\t")[1])
                new_index.append(each_index)
        paired_flag=",".join(temp_flag)#concantate the string
        if paired_flag in pairend_flag_set:#at least this paired end read is mapped some where on the genome
            for_read_info=string.split(string.strip(sam_content[new_index[0]]),"\t")
            rev_read_info=string.split(string.strip(sam_content[new_index[1]]),"\t")
            #quality score for each read
            for_read_qual=for_read_info[4]
            rev_read_qual=rev_read_info[4]
            #mapped chromosome position of each read
            for_read_chromo=for_read_info[2]
            rev_read_chromo=rev_read_info[2]
            record_list=[]
            if for_read_qual > 0 and rev_read_qual > 0 and for_read_chromo==rev_read_chromo:
                query_key=for_read_chromo+";SameChrMapping;Unique"
                if for_read_chromo in ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"]:
                    output_dict[query_key]+=1
    output_dir=sys.argv[3]
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    raw_data_output=open(output_dir+"group"+str(each_group+1)+"_raw_summary.txt",'w')
    for each in output_dict.keys():
        print >>raw_data_output, each, str(output_dict[each])
    raw_data_output.close()
stop = timeit.default_timer()
print stop-start