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
def summarize_mtDNA(data_directory,sample_list):
    for each_group in range(0,40):
        print each_group
        output_dir=data_directory+"all_sample_mtDNA_summary/"
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        output_name="group"+str(each_group+1)+"all_sample_summary.txt"
        output=open(output_dir+"/"+output_name,"w")
        print >>output,"Count","Category","SampleName"
        for each_sample in sample_list:
            print each_sample
            each_sample=string.strip(each_sample)
            input_samfile_name=data_directory+str(each_sample)+"/aligned/group"+str(each_group+1)+"_trimmed_aligned.sam"
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
            else:
                total_group_reads=0
            #count the mapped reads
            correct_noalt_count=0#correctmapping,no alt alignment
            correct_existalt_count=0#correctmapping,with alt alignment
            indistinguishable_count=0#cannot distinguish
            total_count_with_primer=total_group_reads#total number of paired end reads with
            other_count=0
            input_summaryfile_name=data_directory+str(each_sample)+"/summary/group"+str(each_group+1)+"_raw_summary.txt"
            if os.path.isfile(input_summaryfile_name):
                input_summary=open(input_summaryfile_name,"r").readlines()
                for each_line in input_summary:
                    each_line=string.strip(each_line)
                    #remove "(" and ")"
                    each_line=string.replace(each_line,"(","")
                    each_line=string.replace(each_line,")","")
                    each_line=string.replace(each_line,"'","")
                    count_content=string.rsplit(each_line,",",1)
                    mapping_info=string.split(count_content[0],";")
                    count_no=int(count_content[1])
                    alt_info=mapping_info[-1]#must be last element
                    correct_mapping_flag="".join([s for s in mapping_info if "CorrectMapping" in s])
                    if len(correct_mapping_flag)>0:
                        if alt_info=="ExistAltAln":
                            correct_existalt_count=correct_existalt_count+count_no
                        elif alt_info=="NoAltAln":
                            correct_noalt_count=correct_noalt_count+count_no
                    if alt_info=="Indistinguishable":
                        indistinguishable_count=indistinguishable_count+count_no
            print >>output,correct_noalt_count,"CorrectNoAlt",each_sample
            print >>output,correct_existalt_count,"CorrectWithAlt",each_sample
            print >>output,indistinguishable_count,"Indistinguishable",each_sample
            other_count=total_count_with_primer-correct_existalt_count-correct_noalt_count-indistinguishable_count
            print >>output,other_count,"Others",each_sample
        output.close()
import string
import os
exp1_sample=open("C:\Research_data\mtDNA_identification\exp1_list.txt").readlines()
exp2_sample=open("C:\Research_data\mtDNA_identification\exp2_list.txt").readlines()
experiment_group1="exp1_group_result/"
experiment_group2="exp2_group_result/"
data_directory_exp1="C:/Research_data/mtDNA_identification/"+experiment_group1#"C:/Research_data/mtDNA_identification/"
data_directory_exp2="C:/Research_data/mtDNA_identification/"+experiment_group2#"C:/Research_data/mtDNA_identification/"
summarize_mtDNA(data_directory_exp1,exp1_sample)
summarize_mtDNA(data_directory_exp2,exp2_sample)

