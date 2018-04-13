#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=christephen.qin@gmail.com
#SBATCH --mem=50g
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --job-name=mtDNA_remove_primer_and_align_exp1
#SBATCH --output=mtDNAremove_and_align_job_exp1.txt
#SBATCH --partition=general
#read primer data into bash array
module load cutadapt
#pcr_forward_primer=(`cat "pcr_forward_primer.txt"`)
#pcr_reverse_primer=(`cat "pcr_reverse_primer.txt"`)
current_dir=`pwd`
sample_names="NA06224-DTNA_S10
NA06224-ETFB_S4
NA06224-HJ_S22
NA06224-L1PA3_S28
NA06224-mt40n-only_S34
NA18456-ABCC9_S15
NA18456-DTNA_S9
NA18456-ETFB_S3
NA18456-HJ_S21
NA18456-L1PA3_S27
NA18456-mt40n-only_S33
NA18507-ABCC9_S14
NA18507-DTNA_S8
NA18507-ETFB_S2
NA18507-HJ_S20
NA18507-L1PA3_S26
NA18507-mt40n-only_S32"
#NA01359-ABCC9_S13
# NA01359-DTNA_S7
# NA01359-ETFB_S1
# NA01359-HJ_S19
# NA01359-L1PA3_S25
# NA01359-mt40n-only_S31
# NA03330-ABCC9_S17
# NA03330-DTNA_S11
# NA03330-ETFB_S5
# NA03330-HJ_S23
# NA03330-L1PA3_S29
# NA03330-mt40n-only_S35
# NA03576-ABCC9_S18
# NA03576-DTNA_S12
# NA03576-ETFB_S6
# NA03576-HJ_S24
# NA03576-L1PA3_S30
# NA03576-mt40n-only_S36
# NA06224-ABCC9_S16
for sample_name in $sample_names
#for each sample sequenced
do
	#for each group of a primer pair
	genome_file="/home/cp832/qin/bwa.kit/hs38DH.fa"
	input_read1="${current_dir}/exp1_fastq/${sample_name}_L001_R1_001.fastq.gz"
	input_read2="${current_dir}/exp1_fastq/${sample_name}_L001_R2_001.fastq.gz"
	output_trimmed_dir="${current_dir}/exp1_group_result/${sample_name}/trimmed/"
	output_aligned_dir="${current_dir}/exp1_group_result/${sample_name}/aligned/"
	if [ ! -d  "$output_trimmed_dir" ]; then
		mkdir -p "$output_trimmed_dir"
	fi
	if [ ! -d "$output_aligned_dir" ]; then
		mkdir -p "$output_aligned_dir"
	fi	
	for ((i=0;i<=46;i++))
	do
		output_trimmed_file_read1="${current_dir}/exp1_group_result/${sample_name}/trimmed/group$((i+1))_trimmed_read1.fastq"
		output_trimmed_file_read2="${current_dir}/exp1_group_result/${sample_name}/trimmed/group$((i+1))_trimmed_read2.fastq"
		output_aligned_file="${current_dir}/exp1_group_result/${sample_name}/aligned/group$((i+1))_trimmed_aligned.sam"

		#echo $output_aligned_file
		#first remove primers of each group for each paired end reads
		if [ $i -eq 0 ]
		then
			cutadapt -g GGTATGCACTTTTAACAGTCACC -G ACTTGGGTTAATCGTGTGACC --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2
		elif [ $i -eq 1 ]
		then
			cutadapt -g CAGGGTTGGTCAATTTCGT -G TCATAAGGGCTATCGTAGTTTTC --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2		
		elif [ $i -eq 2 ]
		then
			cutadapt -g GTGGCAAGAAATGGGCTAC -G TCTAGTTAATTCATTATGCAGAAGG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2		
		elif [ $i -eq 3 ]
		then
			cutadapt -g AAGCATAATATAGCAAGGACTAACC -G ATGCGGAGGAGAATGTTTTC --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2
		elif [ $i -eq 4 ]
		then
			cutadapt -g AACATATAACTGAACTCCTCACACC -G CCATAGGGTCTTCTCGTCTTG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2
		elif [ $i -eq 5 ]
		then
			cutadapt -g TTTAACCAGTGAAATTGACCTG -G CCTTATTTCTCTTGTCCTTTCGT --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2		
		elif [ $i -eq 6 ]
		then
			cutadapt -g TTTCTATCTACCTTCAAATTCCTC -G GGGTTTTAGGGGCTCTTTG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2		
		elif [ $i -eq 7 ]
		then
			cutadapt -g CTACAACCCTTCGCTGACG -G GTTGGTCTCTGCTAGTGTGGAG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 8 ]
		then
			cutadapt -g GATTACTCCTGCCATCATGAC -G TATCAAAGTAACTCTTTTATCAGACAT --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 9 ]
		then
			cutadapt -g CCATACCCATTACAATCTCCAG -G GATTATGGATGCGGTTGCT --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 10 ]
		then
			cutadapt -g AAACCCTCGTTCCACAGAA -G AATAGTTAAATTAAGAATGGTTATGTT --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 11 ]
		then
			cutadapt -g GATGAATAATAGCAGTTCTACCGT -G GATGAGTGTGGGGAGGAATGG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 12 ]
		then
			cutadapt -g CAACGTAAAAATAAAATGACAGTT -G TGTAAATCTAAAGACAGGGGTTA --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 13 ]
		then
			cutadapt -g ATATGAAAATCACCTCGGAGC -G GGTAGACTGTTCAACCTGTTCC --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 14 ]
		then
			cutadapt -g GCTCGCATCTGCTATAGTGG -G AGTTACAATATGGGAGATTATTCC --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 15 ]
		then
			cutadapt -g ACTATACTACTAACAGACCGCAA -G GGTGAAAAGAAAGATGAATCCTA --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 16 ]
		then
			cutadapt -g AGCATATTTCACCTCCGCTAC -G TATTACTGCTGTTAGAGAAATGAATG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 17 ]
		then
			cutadapt -g ACCCCGATGCATACACCA -G GGACTAGGAAGCAGATAAGGAAA --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 18 ]
		then
			cutadapt -g CGCAAGTAGGTCTACAAGACG -G ATCTGTTTTTAAGCCTAATGTGG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 19 ]
		then
			cutadapt -g CTAATCTTCAACTCCTACATACTTCC -G GCCATACGGTAGTATTTAGTTGG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 20 ]
		then
			cutadapt -g AAGATTAAGAGAACCAACACCTCT -G GTGGCAATAAAAATGATTAAGGA --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 21 ]
		then
			cutadapt -g ACAACACTAAAGGACGAACCTG -G GCATGTGATTGGTGGGTCA --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 22 ]
		then
			cutadapt -g TCATCTTCACAATTCTAATTCTACTG -G CCAGTGCCCTCCTAATTGG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 23 ]
		then
			cutadapt -g GTAACACGAGAAAGCACATACC -G TAGTGAGGAAAGTTGAGCCAAT --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 24 ]
		then
			cutadapt -g TCGAGTCTCCCTTCACCATT -G GGTAAAAGGAGGGCAATTTCT --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 25 ]
		then
			cutadapt -g CGGCTTCGACCCTATATCC -G CTAGTATGGCAATAGGCACAATA --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 26 ]
		then
			cutadapt -g CTACTCTCATAACCCTCAACACC -G ATAGAGAGGTAGAGTTTTTTTCGTG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 27 ]
		then
			cutadapt -g AGCCAACGCCACTTATCC -G AGGGGTAGGCTATGTGTTTTG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 28 ]
		then
			cutadapt -g GCCTCACACTCATTCTCAACC -G TGTAGAGGGAGTATAGGGCTGTG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 29 ]
		then
			cutadapt -g TCTCCTACTTACAGGACTCAACATAC -G ATGCGACAATGGATTTTACATAATGG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 30 ]
		then
			cutadapt -g TTACCACCCTCGTTAACCCT -G GCTGTGTTGGCATCTGCTC --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 31 ]
		then
			cutadapt -g TATATCCTTCTTGCTCATCAGTTG -G GAGTCCTAGTTGACTTGAAGTGG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 32 ]
		then
			cutadapt -g CTTCCACCCCCTAGCAGA -G GTGAGAAGAATTATTCGAGTGCT --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 33 ]
		then
			cutadapt -g AAAGACCACATCATCGAAACC -G AGATAGGGGATTGTGCGGTG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 34 ]
		then
			cutadapt -g CAACATACTCGGATTCTACCCT -G TTGGTGCTGTGGGTGAAAG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 35 ]
		then
			cutadapt -g CCTTCATAAATTATTCAGCTTCCT -G GAGGTCGATGAATGAGTGGTT --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 36 ]
		then
			cutadapt -g CAATGATATGAAAAACCATCGTT -G GGATGGCGGATAGTAAGTTTGT --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 37 ]
		then
			cutadapt -g GCCTATATTACGGATCATTTCTCT -G GGTGTTTAAGGGGTTGGCTA --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 38 ]
		then
			cutadapt -g CTAGGCGACCCAGACAAT -G GCTTCCCCATGAAAGAACA --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2				
		elif [ $i -eq 39 ]
		then
			cutadapt -g CTCCACCATTAGCACCCAAAGC -G TGATTTCACGGAGGATGGTG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2	
		elif [ $i -eq 40 ]
		then
			#EFTB primer
			cutadapt -g TGTGCCATGTTGGTGTGCTG -G ACCCAAAGGATTATAAATCATGCTGCT --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2	
		elif [ $i -eq 41 ]
		then
			#DTNA primer
			cutadapt -g TTCKAGCTTCCCGGCTGCTTT -G CTAGGGAGTGCCAGACAGTGG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2		
		elif [ $i -eq 42 ]
		then
			#HJ primer
			cutadapt -g ACACAGGGAGGGGAACAT -G TGCCATGGTGGTTTGCT --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2
		elif [ $i -eq 43 ]
		then
			#ABCC9 primer
			cutadapt -g TCCCAATGCTATCCCTCCCC -G ATTTGACCCAGCCATCCCAT --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2	
		elif [ $i -eq 44 ]
		then
			#L1PA3 primer
			cutadapt -g GACCCAGCCATCCCATTAC -G CACAACAGTCCCCAGAGTG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2
		elif [ $i -eq 45 ]
		then
			#GAPDH primer
			cutadapt -g CATYRCTRCCACCCAGAA -G GAAWYRAKCYTGACAAAGTGG --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2			
		elif [ $i -eq 46 ]
		then
			#CYCS primer
			cutadapt -g CAAAACACWKKTGAATCTTGC -G TCTCATGACTTTTTYATGTGTACC --pair-filter=any --discard-untrimmed --overlap 10 -o $output_trimmed_file_read1 -p $output_trimmed_file_read2 $input_read1 $input_read2						
		fi
		#then align these paired end reads to the reference genome
		bwa mem -t 10 $genome_file $output_trimmed_file_read1 $output_trimmed_file_read2 > $output_aligned_file
		#then summarize the sam file
	done
done


