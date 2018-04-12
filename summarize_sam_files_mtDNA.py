#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      chris
#
# Created:     29/03/2018
# Copyright:   (c) chris 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import string
import os
from collections import Counter
import sys
input_samfile_name=sys.argv[1]
each_group=int(sys.argv[2])
correct_chr=sys.argv[3]
correct_for_pos=int(sys.argv[4])
correct_rev_pos=int(sys.argv[5])
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
            #further examine which position this paired end is mapped to:
            for_read_info=string.split(string.strip(sam_content[new_index[0]]),"\t")
            rev_read_info=string.split(string.strip(sam_content[new_index[1]]),"\t")
            #quality score for each read
            for_read_qual=for_read_info[4]
            rev_read_qual=rev_read_info[4]
            #mapped chromosome position of each read
            for_read_chromo=for_read_info[2]
            rev_read_chromo=rev_read_info[2]

            #check and record quality score of two mapped paired end reads, a paired end should have a paird of quality scores indicating -log10 Pr(Wrong)
            #a score of "0" means that this read has two mapped positions and their alignment scores are the same and recorded in XA:Z:xxxx field
            #so the summarization stratgy is as follows:
            #1. if NONE of the quality scores are zero: (good mapping)
            # record the quality score and compute genomic position of mapped reads
            #2. if one of the quality score is zero: (one read has multiple mapped position)
            # record the quality score and use the chromosome position of read with non-zero quality score
            # as the optimal alignment on the same chromosome and record possible secondary alignment as well
            #3. if both of the quality score is zero: (two reads have multiple mapped positions on the genome)
            # record the quality score and record all possible alignments on the same chromosome
            #sometimes even mapping quality score is 0, it is still possible that no alternative alignement exists
            ########record list defined here#######
            record_list=[]
            ########record list defined here#######
            #record the quality score of this paired end
            record_list.append(for_read_qual)
            record_list.append(rev_read_qual)
            #check if any alternative alignment exists in mapped reads
            alt_for_read="".join([s for s in for_read_info if "XA:Z:" in s])
            alt_rev_read="".join([s for s in rev_read_info if "XA:Z:" in s])
            #both quality score > 0, suggesting a correct mapped pair end
            if int(for_read_qual)>0 and int(rev_read_qual)>0:
                chrom_for=for_read_chromo
                chrom_for_start="+"+for_read_info[3]
                chrom_for_cigar=for_read_info[5]
                chrom_rev=rev_read_chromo
                chrom_rev_start="-"+rev_read_info[3]
                chrom_rev_cigar=rev_read_info[5]
                record_list.append(chrom_for)
                record_list.append(chrom_for_start)
                record_list.append(chrom_for_cigar)
                record_list.append(chrom_rev)
                record_list.append(chrom_rev_start)
                record_list.append(chrom_rev_cigar)
                if chrom_for==correct_chr and chrom_rev==correct_chr and abs(int(chrom_for_start))>correct_for_pos and abs(int(chrom_rev_start))<correct_rev_pos:
                    record_list.append("CorrectMapping")
                    record_list.append("NoAltAln")
                else:
                    record_list.append("WrongMapping")
                    record_list.append("NoAltAln")
            #Forward read has mapping quality > 0, reverse read might have multiple alignments
            elif int(for_read_qual)>0 and int(rev_read_qual)==0:
                if len(alt_rev_read)>0:#a further check on if the multiple alignment exists in reverse read
                    if for_read_chromo==rev_read_chromo:
                        #Primary alignment
                        chrom_for=for_read_chromo
                        chrom_for_start="+"+for_read_info[3]
                        chrom_for_cigar=for_read_info[5]
                        chrom_rev=rev_read_chromo
                        chrom_rev_start="-"+rev_read_info[3]
                        chrom_rev_cigar=rev_read_info[5]
                        record_list.append(chrom_for)
                        record_list.append(chrom_for_start)
                        record_list.append(chrom_for_cigar)
                        record_list.append(chrom_rev)
                        record_list.append(chrom_rev_start)
                        record_list.append(chrom_rev_cigar)
                        if chrom_for==correct_chr and chrom_rev==correct_chr and abs(int(chrom_for_start))>correct_for_pos and abs(int(chrom_rev_start))<correct_rev_pos:
                            record_list.append("CorrectMapping")
                        else:
                            record_list.append("WrongMapping")
                        #also save alternative alignment of the pair end reads in this case
                        if len(alt_for_read)>0:
                            alt_for_read_chrom=string.split(string.split(alt_for_read,":")[2],";")[0]
                            alt_rev_read_chrom=string.split(string.split(alt_rev_read,":")[2],";")[0]
                            alt_read_combined=alt_for_read_chrom+"/"+alt_rev_read_chrom
                            record_list.append(alt_read_combined)
                            record_list.append("ExistAltAln")
                        else:
                            record_list.append("NoAltAln")
                    else:
                        #e.g read1 chrM 12
                        #read2 chr1 0; XZ:A:chrM
                        chrom_for=for_read_chromo
                        chrom_for_start="+"+for_read_info[3]
                        chrom_for_cigar=for_read_info[5]
                        #because len(alt_rev_read)>0, there must be multiple alignment the in reverse read
                        #we just need to check if alternative alignment could match the chromosome in the forward read
                        alt_rev_info=string.split(string.split(alt_rev_read,":")[2],";")#remove "XA:Z:"
                        alt_rev_aln="".join([s for s in alt_rev_info if chrom_for in s])#find the mapping which has the same chromosome as forward read
                        alt_rev_aln=string.split(alt_rev_aln,",")
                        if len(alt_rev_aln)>1:#a match with the same chromosome of reverse read as in the forward read
                            chrom_rev=alt_rev_aln[0]
                            chrom_rev_start=alt_rev_aln[1]
                            chrom_rev_cigar=alt_rev_aln[2]
                            record_list.append(chrom_for)
                            record_list.append(chrom_for_start)
                            record_list.append(chrom_for_cigar)
                            record_list.append(chrom_rev)
                            record_list.append(chrom_rev_start)
                            record_list.append(chrom_rev_cigar)
                            if chrom_for==correct_chr and chrom_rev==correct_chr and abs(int(chrom_for_start))>correct_for_pos and abs(int(chrom_rev_start))<correct_rev_pos:
                                record_list.append("CorrectMapping")
                            else:
                                record_list.append("WrongMapping")
                        else:
                            record_list.append(chrom_for)
                            record_list.append(chrom_for_start)
                            record_list.append(chrom_for_cigar)
                            record_list.append("NoMatchReverse")
                            record_list.append("NoMapping")
                        if len(alt_for_read)> 0:#forward read with non-zero quality score also has multiple mapping
                            #need to check if these multiple alignments have the same chromosome as in reverse read
                            alt_for_info=string.split(string.split(alt_for_read,":")[2],";")#remove "XA:Z:"
                            alt_for_aln="".join([s for s in alt_for_info if rev_read_chromo in s])#find the mapping which has the same chromosome as reverse read
                            alt_for_aln=string.split(alt_for_aln,",")
                            if len(alt_for_aln)>1:
                                alt_for_read_chrom=alt_for_aln[0]+","+alt_for_aln[1]+","+alt_for_aln[2]
                                alt_rev_read_chrom=rev_read_info[2]+","+"-"+rev_read_info[3]+","+rev_read_info[5]
                                alt_read_combined=alt_for_read_chrom+"/"+alt_rev_read_chrom
                                record_list.append(alt_read_combined)
                                record_list.append("ExistAltAln")
                            else:
                                alt_rev_read_chrom=rev_read_info[2]+","+"-"+rev_read_info[3]+","+rev_read_info[5]
                                alt_read_combined="NoAltMatchForward"+"/"+alt_rev_read_chrom
                                record_list.append(alt_read_combined)
                                record_list.append("NoAltAln")
                        else:
                            record_list.append("NoAltAln")
                else:#although quality score of reverse read is 0, no alternative alignment exists
                        chrom_for=for_read_chromo
                        chrom_for_start="+"+for_read_info[3]
                        chrom_for_cigar=for_read_info[5]
                        chrom_rev=rev_read_chromo
                        chrom_rev_start="-"+rev_read_info[3]
                        chrom_rev_cigar=rev_read_info[5]
                        record_list.append(chrom_for)
                        record_list.append(chrom_for_start)
                        record_list.append(chrom_for_cigar)
                        record_list.append(chrom_rev)
                        record_list.append(chrom_rev_start)
                        record_list.append(chrom_rev_cigar)
                        if chrom_for==correct_chr and chrom_rev==correct_chr and abs(int(chrom_for_start))>correct_for_pos and abs(int(chrom_rev_start))<correct_rev_pos:
                            record_list.append("CorrectMapping")
                        else:
                            record_list.append("WrongMapping")
                        record_list.append("NoAltAln")
            #Reverse read has mapping quality >0, forward read might have multiple alignments
            elif int(for_read_qual)==0 and int(rev_read_qual)>0:
                if len(alt_for_read)>0:
                    if for_read_chromo==rev_read_chromo:
                        #Primary alignment
                        chrom_for=for_read_chromo
                        chrom_for_start="+"+for_read_info[3]
                        chrom_for_cigar=for_read_info[5]
                        chrom_rev=rev_read_chromo
                        chrom_rev_start="-"+rev_read_info[3]
                        chrom_rev_cigar=rev_read_info[5]
                        record_list.append(chrom_for)
                        record_list.append(chrom_for_start)
                        record_list.append(chrom_for_cigar)
                        record_list.append(chrom_rev)
                        record_list.append(chrom_rev_start)
                        record_list.append(chrom_rev_cigar)
                        #also save alternative alignment of the pair end reads in this case
                        if chrom_for==correct_chr and chrom_rev==correct_chr and abs(int(chrom_for_start))>correct_for_pos and abs(int(chrom_rev_start))<correct_rev_pos:
                            record_list.append("CorrectMapping")
                        else:
                            record_list.append("WrongMapping")
                        if len(alt_rev_read)>0:
                            alt_for_read_chrom=string.split(string.split(alt_for_read,":")[2],";")[0]
                            alt_rev_read_chrom=string.split(string.split(alt_rev_read,":")[2],";")[0]
                            alt_read_combined=alt_for_read_chrom+"/"+alt_rev_read_chrom
                            record_list.append(alt_read_combined)
                            record_list.append("ExistAltAln")
                        else:
                            record_list.append("NoAltAln")
                    else:
                        #different chromosome mapping in the primary alignment
                        #e.g read1 chr1 0, XA:Z:
                        #read2 chrM 12
                        chrom_rev=rev_read_chromo
                        chrom_rev_start="-"+rev_read_info[3]
                        chrom_rev_cigar=rev_read_info[5]
                        if len(alt_for_read)>0:#sometimes even a quality score of 0 does not have multiple alignment
                            alt_for_info=string.split(string.split(alt_for_read,":")[2],";")#remove "XA:Z:"
                            alt_for_aln="".join([s for s in alt_for_info if rev_read_chromo in s])#find the mapping which has the same chromosome as forward read
                            alt_for_aln=string.split(alt_for_aln,",")
                            #split the information and save it
                            if len(alt_for_aln)>1:#there is a match in the multiple alignments
                                chrom_for=alt_for_aln[0]
                                chrom_for_start=alt_for_aln[1]
                                chrom_for_cigar=alt_for_aln[2]
                                record_list.append(chrom_for)
                                record_list.append(chrom_for_start)
                                record_list.append(chrom_for_cigar)
                                record_list.append(chrom_rev)
                                record_list.append(chrom_rev_start)
                                record_list.append(chrom_rev_cigar)
                                if chrom_for==correct_chr and chrom_rev==correct_chr and abs(int(chrom_for_start))>correct_for_pos and abs(int(chrom_rev_start))<correct_rev_pos:
                                    record_list.append("CorrectMapping")
                                else:
                                    record_list.append("WrongMapping")
                            else:
                                record_list.append("NoMatchForward")
                                record_list.append(chrom_rev)
                                record_list.append(chrom_rev_start)
                                record_list.append(chrom_rev_cigar)
                                record_list.append("NoMapping")
                            if len(alt_rev_read)> 0:#reverse read of non-zero quality score also has multiple mapping
                                #look for alternative alignment with same chromosome mapping
                                alt_rev_info=string.split(string.split(alt_rev_read,":")[2],";")#remove "XA:Z:"
                                alt_rev_aln="".join([s for s in alt_rev_info if for_read_chromo in s])#find the mapping which has the same chromosome as reverse read
                                alt_rev_aln=string.split(alt_rev_aln,",")
                                if len(alt_rev_aln)>1:#find a match which has the same chromosome as in forward read
                                    alt_rev_read_chrom=alt_rev_aln[0]+","+alt_rev_aln[1]+","+alt_rev_aln[2]
                                    alt_for_read_chrom=for_read_info[2]+","+"+"+for_read_info[3]+","+for_read_info[5]
                                    alt_read_combined=alt_for_read_chrom+"/"+alt_rev_read_chrom
                                    record_list.append(alt_read_combined)
                                    record_list.append("ExistAltAln")
                                else:
                                    alt_for_read_chrom=for_read_info[2]+","+"+"+for_read_info[3]+","+for_read_info[5]
                                    alt_read_combined=alt_for_read_chrom+"/"+"NoAltMatchReverse"
                                    record_list.append(alt_read_combined)
                                    record_list.append("NoAltAln")
                            else:
                                record_list.append("NoAltAln")
                else:#although quality score of reverse read is 0, no alternative alignment exists
                        chrom_for=for_read_chromo
                        chrom_for_start="+"+for_read_info[3]
                        chrom_for_cigar=for_read_info[5]
                        chrom_rev=rev_read_chromo
                        chrom_rev_start="-"+rev_read_info[3]
                        chrom_rev_cigar=rev_read_info[5]
                        record_list.append(chrom_for)
                        record_list.append(chrom_for_start)
                        record_list.append(chrom_for_cigar)
                        record_list.append(chrom_rev)
                        record_list.append(chrom_rev_start)
                        record_list.append(chrom_rev_cigar)
                        if chrom_for==correct_chr and chrom_rev==correct_chr and abs(int(chrom_for_start))>correct_for_pos and abs(int(chrom_rev_start))<correct_rev_pos:
                            record_list.append("CorrectMapping")
                        else:
                            record_list.append("WrongMapping")
                        record_list.append("NoAltAln")
            else:
                #this situation means there are alternative alignments in both reads
                #and cannot be distinguished
                if len(alt_for_read)>0 and len(alt_rev_read)>0:
                    if for_read_chromo==rev_read_chromo:
                        #Primary alignment
                        chrom_for=for_read_chromo
                        chrom_for_start="+"+for_read_info[3]
                        chrom_for_cigar=for_read_info[5]
                        chrom_rev=rev_read_chromo
                        chrom_rev_start="-"+rev_read_info[3]
                        chrom_rev_cigar=rev_read_info[5]
                        record_list.append(chrom_for)
                        record_list.append(chrom_for_start)
                        record_list.append(chrom_for_cigar)
                        record_list.append(chrom_rev)
                        record_list.append(chrom_rev_start)
                        record_list.append(chrom_rev_cigar)
                        alt_for_read_chrom=string.split(string.split(alt_for_read,":")[2],";")[0]
                        alt_rev_read_chrom=string.split(string.split(alt_rev_read,":")[2],";")[0]
                        alt_read_combined=alt_for_read_chrom+"/"+alt_rev_read_chrom
                        record_list.append(alt_read_combined)
                        record_list.append("Indistinguishable")
                    else:
                        chrom_for=for_read_chromo
                        chrom_for_start="+"+for_read_info[3]
                        chrom_for_cigar=for_read_info[5]
                        alt_rev_info=string.split(string.split(alt_rev_read,":")[2],";")#remove "XA:Z:"
                        alt_rev_aln="".join([s for s in alt_rev_info if chrom_for in s])#find the mapping which has the same chromosome as forward read
                        alt_rev_aln=string.split(alt_rev_aln,",")
                        #split the information and save it
                        if len(alt_rev_aln)>1:#a match with the same chromosome of reverse read as in the forward read
                            chrom_rev=alt_rev_aln[0]
                            chrom_rev_start=alt_rev_aln[1]
                            chrom_rev_cigar=alt_rev_aln[2]
                            record_list.append(chrom_for)
                            record_list.append(chrom_for_start)
                            record_list.append(chrom_for_cigar)
                            record_list.append(chrom_rev)
                            record_list.append(chrom_rev_start)
                            record_list.append(chrom_rev_cigar)
                        else:
                            record_list.append(chrom_for)
                            record_list.append(chrom_for_start)
                            record_list.append(chrom_for_cigar)
                            record_list.append("NoMatchReverse")
                        record_list.append("or")
                        chrom_rev_alt=rev_read_chromo
                        chrom_rev_start_alt=rev_read_info[2]
                        chrom_rev_cigar_alt=rev_read_info[5]
                        alt_for_info=string.split(string.split(alt_for_read,":")[2],";")#remove "XA:Z:"
                        alt_for_aln="".join([s for s in alt_for_info if rev_read_chromo in s])#find the mapping which has the same chromosome as forward read
                        alt_for_aln=string.split(alt_for_aln,",")
                        if len(alt_for_aln)>1:
                            chrom_for_alt=alt_for_aln[0]
                            chrom_for_start_alt=alt_for_aln[1]
                            chrom_for_cigar_alt=alt_for_aln[2]
                            record_list.append(chrom_for)
                            record_list.append(chrom_for_start)
                            record_list.append(chrom_for_cigar)
                            record_list.append(chrom_rev)
                            record_list.append(chrom_rev_start)
                            record_list.append(chrom_rev_cigar)
                        else:
                            record_list.append("NoMatchForward")
                            record_list.append(chrom_rev_alt)
                            record_list.append(chrom_rev_start_alt)
                            record_list.append(chrom_rev_cigar_alt)
                        record_list.append("Indistinguishable")
                else:#at least one read has no multiple alignment sequence
                    chrom_for=for_read_chromo
                    chrom_for_start="+"+for_read_info[3]
                    chrom_for_cigar=for_read_info[5]
                    chrom_rev=rev_read_chromo
                    chrom_rev_start="-"+rev_read_info[3]
                    chrom_rev_cigar=rev_read_info[5]
                    record_list.append(chrom_for)
                    record_list.append(chrom_for_start)
                    record_list.append(chrom_for_cigar)
                    record_list.append(chrom_rev)
                    record_list.append(chrom_rev_start)
                    record_list.append(chrom_rev_cigar)
                    if chrom_for==correct_chr and chrom_rev==correct_chr and abs(int(chrom_for_start))>correct_for_pos and abs(int(chrom_rev_start))<correct_rev_pos:
                        record_list.append("CorrectMapping")
                    else:
                        record_list.append("WrongMapping")
                    record_list.append("NoAltAln")
            #finally save record_list into the sam summary list
            sam_summary_all.append(";".join(record_list))
    sam_summary_count=Counter(sam_summary_all)
    output_dir=sys.argv[6]
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    raw_data_output=open(output_dir+"group"+str(each_group+1)+"_raw_summary.txt",'w')
    for each in sam_summary_count.most_common():
        print >>raw_data_output, each
    raw_data_output.close()
else:
    output_dir=sys.argv[6]
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    raw_data_output=open(output_dir+"group"+str(each_group+1)+"_raw_summary.txt",'w')
    print >>raw_data_output, "No reads found."
    raw_data_output.close()
##        raw_content=[]
##        summary_content=[]
##
##        read_head_summary=Counter(all_read_head)
##        if read_head_summary.most_common()[0][1]>2:
##            print read_head_summary.most_common()[0:5]

##            raw_content.append(str(read1)+str(read2))
##            chr_read1=read1[1]
##            chr_read2=read2[1]
##            sam_flag_read1=read1[0]
##            sam_flag_read2=read2[0]
##            if chr_read1=="chrM" and chr_read2=="chrM":
##                start_position=read1[2]
##                end_position=str(int(read1[2])+int(read1[7])-1)
##                combined=sam_flag_read1+","+sam_flag_read2+","+chr_read1+","+start_position+","+end_position
##            else:
##                combined=sam_flag_read1+","+sam_flag_read2+","+chr_read1+","+chr_read2
##            summary_content.append(combined)
##        raw_summary=Counter(raw_content)
##        detail_summary=Counter(summary_content)
##        raw_output=open("C:/Research_data/mtDNA_identification/exp1_group_result/"+str(each_sample)+"/summary/group"+str(i+1)+"raw_summary.txt","w")
##        detail_output=open("C:/Research_data/mtDNA_identification/exp1_group_result/"+str(each_sample)+"/summary/group"+str(i+1)+"detail_summary.txt","w")
##        for each in raw_summary.most_common():
##            print >>raw_output,str(each[0])+","+str(each[1])
##        raw_output.close()
##        for each in detail_summary.most_common():
##            print >>detail_output, str(each[0])+","+str(each[1])
##        detail_output.close()
##for each_sample in exp2_sample:
##    for i in range(0,43):
##        input_samfile_name="C:/Research_data/mtDNA_identification/exp2_group_result/"+str(each_sample)+"/aligned/group"+str(i+1)+"_trimmed_aligned.sam"
##        input_samfile=open(input_samfile_name,"r").readlines()
##        sam_content=[]
##        for each_line in input_samfile:
##            if each_line[0]!="@":
##                sam_content.append(each_line)
##        raw_content=[]
##        summary_content=[]
##        for j in range(0,len(sam_content),2):
##            read1=string.split(sam_content[j],"\t")[1:9]
##            read2=string.split(sam_content[j+1],"\t")[1:9]
##            raw_content.append(str(read1)+str(read2))
##            chr_read1=read1[1]
##            chr_read2=read2[1]
##            sam_flag_read1=read1[0]
##            sam_flag_read2=read2[0]
##            if chr_read1=="chrM" and chr_read2=="chrM":
##                start_position=read1[2]
##                end_position=str(int(read1[2])+int(read1[7]))
##                combined=sam_flag_read1+","+sam_flag_read2+","+chr_read1+","+start_position+","+end_position
##            else:
##                combined=sam_flag_read1+","+sam_flag_read2+","+chr_read1+","+chr_read2
##            summary_content.append(combined)
##        raw_summary=Counter(raw_content)
##        detail_summary=Counter(summary_content)
##        raw_output=open("C:/Research_data/mtDNA_identification/exp2_group_result/"+str(each_sample)+"/group"+str(i+1)+"/raw_summary.txt","w")
##        detail_output=open("C:/Research_data/mtDNA_identification/exp2_group_result/"+str(each_sample)+"/group"+str(i+1)+"/detail_summary.txt","w")
##        for each in raw_summary.most_common():
##            print >>raw_output,str(each[0])+","+str(each[1])
##        raw_output.close()
##        for each in detail_summary.most_common():
##            print >>detail_output, str(each[0])+","+str(each[1])
##        detail_output.close()
