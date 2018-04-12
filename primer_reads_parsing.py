#-------------------------------------------------------------------------------
# Name:        primer reads parsing
# Purpose: This script aims at parsing reads derived with use of 40 pairs of PCR primers and write them into separate fastq files for further mapping
#
# Author:      chris
#
# Created:     21/03/2018
# Copyright:   (c) chris 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import string
from Bio.Seq import Seq
def check_primer_group(input_read1,input_read2,pcr_primer_seq_db):
    forward_check=0
    reverse_check=0
    group_check=0
    perfect_check=0
    for each_primer in range(0,len(pcr_primer_seq_db)):
        [pcr_for_seq,pcr_rev_seq]=string.split(pcr_primer_seq_db[each_primer],sep=",")
        if input_read1[0:len(pcr_for_seq)]==pcr_for_seq:
            forward_check=1
        if input_read2[0:len(pcr_rev_seq)]==pcr_rev_seq:
            reverse_check=1
        if forward_check==1 or reverse_check==1:
            if forward_check==1 and reverse_check==1:
                perfect_check=1
            group_check=each_primer+1
            break
    output=[forward_check,reverse_check,group_check,perfect_check]
    return output
#read fastq file paired end reads into one list
pcr_primer_sequence_file=open("C:\Research_data\mtDNA_identification\primer_40_seq.csv","r").readlines()
#first convert the PCR primer pair to the sequence easier to compare
pcr_primer_sequence_database=list()
for each_primer in range(0,len(pcr_primer_sequence_file)):
    primer_info=string.split(string.strip(pcr_primer_sequence_file[each_primer]),sep=",")#split PCR primer line into different fields
    primer_forward=Seq(primer_info[1])#forward PCR primer sequence
    primer_reverse=Seq(primer_info[3])#reverse PCR primer sequence
    primer_reverse_comp=primer_reverse.reverse_complement()#then reverse the sequence and get the complement
    primer_all=str(primer_forward)+","+str(primer_reverse_comp)
    pcr_primer_sequence_database.append(primer_all)
pcr_primer_output=open("C:\\Research_data\\mtDNA_identification\\pcr_primer_sequence.txt",'w')
for each in pcr_primer_sequence_database:
    print >>pcr_primer_output,each
pcr_primer_output.close()
#get each paired end reads, find which pcr primer it belongs to and write to the corresponding files
#sample_name="NA18456-mt40n-only_S33"
for sample_name in ["NA01359-mt40n-only_S31","NA03330-mt40n-only_S35","NA06224-mt40n-only_S34","NA18507-mt40n-only_S32"]:
    print sample_name
    input_read1_file="C:\\Research_data\\mtDNA_identification\\exp1\\mt40n_only\\"+sample_name+"_L001_R1_001.fastq"
    input_read2_file="C:\\Research_data\\mtDNA_identification\\exp1\\mt40n_only\\"+sample_name+"_L001_R2_001.fastq"
    input_read1=open(input_read1_file,"r").readlines()
    input_read2=open(input_read2_file,"r").readlines()
    for each_pe in range(0,len(input_read1),4):
        if each_pe/len(input_read1)>0.2:
            print "20%"
        elif each_pe/len(input_read1)>0.5:
            print "50%"

        forward_read=string.strip(input_read1[each_pe+1])#forward read 1
        reverse_read=string.strip(input_read2[each_pe+1])#reverse read 2
        [for_check,reverse_check,group_id,perfect]=check_primer_group(forward_read,reverse_read,pcr_primer_sequence_database)
        if perfect==1:
            output_file_name1="C:\\Research_data\\mtDNA_identification\\group_result\\"+sample_name+"\\"+str(group_id)+"p_read1.fastq"
            output_file_name2="C:\\Research_data\\mtDNA_identification\\group_result\\"+sample_name+"\\"+str(group_id)+"p_read2.fastq"
        else:
            output_file_name1="C:\\Research_data\\mtDNA_identification\\group_result\\"+sample_name+"\\"+str(group_id)+"_read1.fastq"
            output_file_name2="C:\\Research_data\\mtDNA_identification\\group_result\\"+sample_name+"\\"+str(group_id)+"_read2.fastq"
        output1=open(output_file_name1,'a')
        output2=open(output_file_name2,'a')
        print >>output1,string.strip(input_read1[each_pe])+'\t'+str(for_check)+','+str(reverse_check)+','+str(group_id)
        print >>output1,input_read1[each_pe+1],
        print >>output1,input_read1[each_pe+2],
        print >>output1,input_read1[each_pe+3],
        print >>output2,string.strip(input_read2[each_pe])+'\t'+str(for_check)+','+str(reverse_check)+','+str(group_id)
        print >>output2,input_read2[each_pe+1],
        print >>output2,input_read2[each_pe+2],
        print >>output2,input_read2[each_pe+3],
        output1.close()
        output2.close()
##sample_idx=2264245
##sample_pe_read1=string.strip(input_read1[sample_idx])
##sample_pe_read2=string.strip(input_read2[sample_idx])
##sample_group_id=check_primer_group(sample_pe_read1,sample_pe_read2,pcr_primer_sequence_database)
##print sample_group_id