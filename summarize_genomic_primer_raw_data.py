#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      chris
#
# Created:     09/04/2018
# Copyright:   (c) chris 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------
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
def summarize_genomicDNA(data_directory,sample_list,each_group):
    #group 41: ETFB
    #group 42: DTNA
    #group 43: HJ
    #group 44: ABCC9
    #group 45: L1PA3
    print each_group
    output_dir=data_directory+"all_sample_genomicDNA_summary/"
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    output_name="group"+str(each_group+1)+"all_sample_summary.txt"
    output=open(output_dir+"/"+output_name,"w")
    print >>output,"Count","Chr","SampleName"
    for each_sample in sample_list:
        each_sample=string.strip(each_sample)
        print each_sample
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
        chr1_count=0
        chr2_count=0
        chr3_count=0
        chr4_count=0
        chr5_count=0
        chr6_count=0
        chr7_count=0
        chr8_count=0
        chr9_count=0
        chr10_count=0
        chr11_count=0
        chr12_count=0
        chr13_count=0
        chr14_count=0
        chr15_count=0
        chr16_count=0
        chr17_count=0
        chr18_count=0
        chr19_count=0
        chr20_count=0
        chr21_count=0
        chr22_count=0
        chrX_count=0
        chrY_count=0
        chrM_count=0
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
                #mapping_info=string.split(count_content[0],";")
                count_no=int(count_content[1])
                count_content=count_content[0]
                if count_content=="chr1;SameChrMapping;Unique":
                    chr1_count=count_no
                elif count_content=="chr2;SameChrMapping;Unique":
                    chr2_count=count_no
                elif count_content=="chr3;SameChrMapping;Unique":
                    chr3_count=count_no
                elif count_content=="chr4;SameChrMapping;Unique":
                    chr4_count=count_no
                elif count_content=="chr5;SameChrMapping;Unique":
                    chr5_count=count_no
                elif count_content=="chr6;SameChrMapping;Unique":
                    chr6_count=count_no
                elif count_content=="chr7;SameChrMapping;Unique":
                    chr7_count=count_no
                elif count_content=="chr8;SameChrMapping;Unique":
                    chr8_count=count_no
                elif count_content=="chr9;SameChrMapping;Unique":
                    chr9_count=count_no
                elif count_content=="chr10;SameChrMapping;Unique":
                    chr10_count=count_no
                elif count_content=="chr11;SameChrMapping;Unique":
                    chr11_count=count_no
                elif count_content=="chr12;SameChrMapping;Unique":
                    chr12_count=count_no
                elif count_content=="chr13;SameChrMapping;Unique":
                    chr13_count=count_no
                elif count_content=="chr14;SameChrMapping;Unique":
                    chr14_count=count_no
                elif count_content=="chr15;SameChrMapping;Unique":
                    chr15_count=count_no
                elif count_content=="chr16;SameChrMapping;Unique":
                    chr16_count=count_no
                elif count_content=="chr17;SameChrMapping;Unique":
                    chr17_count=count_no
                elif count_content=="chr18;SameChrMapping;Unique":
                    chr18_count=count_no
                elif count_content=="chr19;SameChrMapping;Unique":
                    chr19_count=count_no
                elif count_content=="chr20;SameChrMapping;Unique":
                    chr20_count=count_no
                elif count_content=="chr21;SameChrMapping;Unique":
                    chr21_count=count_no
                elif count_content=="chr22;SameChrMapping;Unique":
                    chr22_count=count_no
                elif count_content=="chr3;SameChrMapping;Unique":
                    chr3_count=count_no
                elif count_content=="chr3;SameChrMapping;Unique":
                    chr3_count=count_no
                elif count_content=="chrX;SameChrMapping;Unique":
                    chrX_count=count_no
                elif count_content=="chrY;SameChrMapping;Unique":
                    chrY_count=count_no
                elif count_content=="chrM;SameChrMapping;Unique":
                    chrM_count=count_no
            other_count=total_group_reads-chr1_count-chr2_count-chr3_count-chr4_count-chr5_count-chr6_count-chr7_count-chr8_count-chr9_count-chr10_count-chr11_count-chr12_count-chr13_count-chr14_count-chr15_count-chr16_count-chr17_count-chr18_count-chr19_count-chr20_count-chr21_count-chr22_count-chrX_count-chrY_count-chrM_count
            print >>output, chr1_count,"chr1",each_sample
            print >>output, chr2_count,"chr2",each_sample
            print >>output, chr3_count,"chr3",each_sample
            print >>output, chr4_count,"chr4",each_sample
            print >>output, chr5_count,"chr5",each_sample
            print >>output, chr6_count,"chr6",each_sample
            print >>output, chr7_count,"chr7",each_sample
            print >>output, chr8_count,"chr8",each_sample
            print >>output, chr9_count,"chr9",each_sample
            print >>output, chr10_count,"chr10",each_sample
            print >>output, chr11_count,"chr11",each_sample
            print >>output, chr12_count,"chr12",each_sample
            print >>output, chr13_count,"chr13",each_sample
            print >>output, chr14_count,"chr14",each_sample
            print >>output, chr15_count,"chr15",each_sample
            print >>output, chr16_count,"chr16",each_sample
            print >>output, chr17_count,"chr17",each_sample
            print >>output, chr18_count,"chr18",each_sample
            print >>output, chr19_count,"chr19",each_sample
            print >>output, chr20_count,"chr20",each_sample
            print >>output, chr21_count,"chr21",each_sample
            print >>output, chr22_count,"chr22",each_sample
            print >>output, chrX_count,"chrX",each_sample
            print >>output, chrY_count,"chrY",each_sample
            print >>output, chrM_count,"chrM",each_sample
            print >>output, other_count,"other",each_sample
        else:
            print >>output, chr1_count,"chr1",each_sample
            print >>output, chr2_count,"chr2",each_sample
            print >>output, chr3_count,"chr3",each_sample
            print >>output, chr4_count,"chr4",each_sample
            print >>output, chr5_count,"chr5",each_sample
            print >>output, chr6_count,"chr6",each_sample
            print >>output, chr7_count,"chr7",each_sample
            print >>output, chr8_count,"chr8",each_sample
            print >>output, chr9_count,"chr9",each_sample
            print >>output, chr10_count,"chr10",each_sample
            print >>output, chr11_count,"chr11",each_sample
            print >>output, chr12_count,"chr12",each_sample
            print >>output, chr13_count,"chr13",each_sample
            print >>output, chr14_count,"chr14",each_sample
            print >>output, chr15_count,"chr15",each_sample
            print >>output, chr16_count,"chr16",each_sample
            print >>output, chr17_count,"chr17",each_sample
            print >>output, chr18_count,"chr18",each_sample
            print >>output, chr19_count,"chr19",each_sample
            print >>output, chr20_count,"chr20",each_sample
            print >>output, chr21_count,"chr21",each_sample
            print >>output, chr22_count,"chr22",each_sample
            print >>output, chrX_count,"chrX",each_sample
            print >>output, chrY_count,"chrY",each_sample
            print >>output, chrM_count,"chrM",each_sample
            print >>output, other_count,"other",each_sample
    output.close()
import string
import os
exp1_sample=open("C:\Research_data\mtDNA_identification\exp1_list.txt").readlines()
exp2_sample=open("C:\Research_data\mtDNA_identification\exp2_list.txt").readlines()
experiment_group1="exp1_group_result/"
experiment_group2="exp2_group_result/"
data_directory_exp1="C:/Research_data/mtDNA_identification/"+experiment_group1#"C:/Research_data/mtDNA_identification/"
data_directory_exp2="C:/Research_data/mtDNA_identification/"+experiment_group2#"C:/Research_data/mtDNA_identification/"
#group 41: ETFB
#group 42: DTNA
#group 43: HJ
#group 44: ABCC9
#group 45: L1PA3
genomic_primer_list_exp1=["ETFB","DTNA","HJ","ABCC9","L1PA3"]
genomic_primer_list_exp2=["ETFB","DTNA","HJ"]
#exp1/exp2
##for each_primer in genomic_primer_list_exp1:
##    all_samples=[s for s in exp1_sample if each_primer in s]#match this primer
##    if each_primer=="ETFB":
##        group=40
##    elif each_primer=="DTNA":
##        group=41
##    elif each_primer=="HJ":
##        group=42
##    elif each_primer=="ABCC9":
##        group=43
##    elif each_primer=="L1PA3":
##        group=44
##    summarize_genomicDNA(data_directory_exp2,all_samples,group)
###exp2
for each_primer in genomic_primer_list_exp2:
    all_samples=[s for s in exp2_sample if each_primer in s]#match this primer
    all_samples.append("NA01359-mt40n-all_S29")
    all_samples.append("NA03330-mt40n-all_S33")
    all_samples.append("NA03576-mt40n-all_S34")
    all_samples.append("NA06224-mt40n-all_S32")
    all_samples.append("NA12156-mt40n-all_S35")
    all_samples.append("NA18456-mt40n-all_S31")
    all_samples.append("NA18507-mt40n-all_S30")
    if each_primer=="ETFB":
        group=40
    elif each_primer=="DTNA":
        group=41
    elif each_primer=="HJ":
        group=42
    summarize_genomicDNA(data_directory_exp2,all_samples,group)
#CYCS-GAPDH primer in exp2
genomic_primer_gapdh_cycs=["GAPDH","CYCS"]
for each_primer in genomic_primer_gapdh_cycs:
    all_samples=[s for s in exp2_sample if each_primer in s]#match this primer
    all_samples.append("NA01359-mt40n-all_S29")
    all_samples.append("NA03330-mt40n-all_S33")
    all_samples.append("NA03576-mt40n-all_S34")
    all_samples.append("NA06224-mt40n-all_S32")
    all_samples.append("NA12156-mt40n-all_S35")
    all_samples.append("NA18456-mt40n-all_S31")
    all_samples.append("NA18507-mt40n-all_S30")
    if each_primer=="GAPDH":
        group=45
    elif each_primer=="CYCS":
        group=46
    summarize_genomicDNA(data_directory_exp2,all_samples,group)