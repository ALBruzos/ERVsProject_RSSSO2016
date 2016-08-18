#!/bin/bash
#######################################################################
#            PROJECT: Endogenous retroviruses in Bos Taurus           #
#                       Subtask: 8.1 ILLUMINA DATA v2                 #
#             Research Summer School in Statistical Omics             #
#                     author: Alicia L. Bruzos                        #
#                      date: 16th August 2016                         #
#######################################################################
#Repeat the script 6.2 for the illumina data, so add a conversion from bam files to fastq files.

#set -e  #Stop the script if it has an error

##Step 1: get all the information that we will need: Positions lists of ERV2, BTLTR, fastq files, BWA, samtools, reference genomes.
BWA=/common/WORK/SCHOOL2016/Data/ERVs/bwa/bwa
samtools=/common/WORK/SCHOOL2016/Data/ERVs/samtools-1.3.1/samtools

ERV2ref=/common/WORK/SCHOOL2016/Data/ERVs/refs/ERV2_Full.fa
BTLref=/common/WORK/SCHOOL2016/Data/ERVs/refs/BTLTR1_ref_v1.fa

BBBreads=/common/WORK/SCHOOL2016/Data/ERVs/BBB/

##Step 2: create a list of the .bam files containing the illumina reads 
ls ${BBBreads}*.bam > bamfiles_BBB.txt
echo "bamfiles_BBB.txt has been created. Look in this path: /common/WORK/SCHOOL2016/students/alicia00/Step7/"
BBBfiles=/common/WORK/SCHOOL2016/students/alicia00/Step7/bamfiles_BBB.txt

		#Results bamfiles_BBB.txt:
		# /common/WORK/SCHOOL2016/Data/ERVs/BBB/chr10_102477536.discordant.bam
		# /common/WORK/SCHOOL2016/Data/ERVs/BBB/chr10_33763278.discordant.bam
		# /common/WORK/SCHOOL2016/Data/ERVs/BBB/chr10_34445397.discordant.bam
		# /common/WORK/SCHOOL2016/Data/ERVs/BBB/chr10_42546471.discordant.bam

##Step 3: Create my new positions lists (only start position!!) for ERV2 and LTR from BBBdataset.
#done manually with excell
ST_ERV2_BBB=/common/WORK/SCHOOL2016/Data/ERVs/examples/ERV2_starts.txt
ST_BTL_BBB=/common/WORK/SCHOOL2016/Data/ERVs/examples/BTLTR1_starts.txt
#ST_BTL_BBB=/common/WORK/SCHOOL2016/students/alicia00/Step7/ST_BTL_BBB.txt
#Create bash arrays for the position lists of ERV2 and BTLTR
ERV2pl=(`cat $ST_ERV2_BBB`) #getting all the lines of that file into a bash array
BTLpl=(`cat $ST_BTL_BBB`)

##Step 4: Create two for loops to apply the following steps to ERV2 and to BTL...
#For loop for ERV2
for itemA in ${ERV2pl[@]}  #for each item of that array do the following
	do
		echo "ERV2 position tested: ${itemA}"
		FILENAME_path=`grep ${itemA} ${BBBfiles}`
		##/common/WORK/SCHOOL2016/Data/ERVs/BBB/chr10_102477536.discordant.bam
		echo "My input is ${FILENAME_path}"
		FILENAME_file=`echo ${FILENAME_path} | awk '{n=split($0,end_arr,"/"); print end_arr[n]}'`
		${samtools} sort -n  ${FILENAME_path} | ${samtools} fastq -1 ${FILENAME_file}.r1.fq -2 ${FILENAME_file}.r2.fq -
		echo "fastq files have been generated for ${itemA}"
		
		if [ "${FILENAME_file}" == "" ]
			then
				echo "There are not fastq files corresponding with the position ${itemA}"
			
		else if [ "${FILENAME_file}" != "" ]
			then
				#Running BWAMEM to align in that file and getting the sam output file (changing for illumina: https://github.com/lh3/bwa)
				${BWA} mem -t 6 ${ERV2ref} ${FILENAME_file}.r1.fq ${FILENAME_file}.r2.fq > ${FILENAME_file}.aligned.sam
				echo "BWA has run for ${FILENAME_file}"
				
				#Running samtools to convert sam file to bam file
				${samtools} view -@ 6 -b ${FILENAME_file}.aligned.sam > ${FILENAME_file}.aligned.bam
				
				#Running samtools to sort the files
				${samtools} sort -@ 6 -o ${FILENAME_file}.aligned.sorted.bam ${FILENAME_file}.aligned.bam

				#Running samtools to index the files
				${samtools} index ${FILENAME_file}.aligned.sorted.bam
				echo "Samtools (sam-bam conversion, sort and index) has run for ${FILENAME_file}"
			fi
			fi
	done

echo "END OF ANALYSIS OF ERV2"


#For loop for BTL
for itemB in ${BTLpl[@]}  #for each item of that array do the following
	do
		echo "BTLTR position tested: ${itemB}"
		FILENAME_path=`grep ${itemB} ${BBBfiles}`
		##/common/WORK/SCHOOL2016/Data/ERVs/BBB/chr10_102477536.discordant.bam
		echo "My input is ${FILENAME_path}"
		FILENAME_file=`echo ${FILENAME_path} | awk '{n=split($0,end_arr,"/"); print end_arr[n]}'`
		${samtools} sort -n  ${FILENAME_path} | ${samtools} fastq -1 ${FILENAME_file}.r1.fq -2 ${FILENAME_file}.r2.fq -
		echo "fastq files have been generated for ${itemB}"
		
		if [ "${FILENAME_file}" == "" ]
			then
				echo "There are not fastq files corresponding with the position ${itemA}"
			
		else if [ "${FILENAME_file}" != "" ]
			then
				#Running BWAMEM to align in that file and getting the sam output file (changing for illumina: https://github.com/lh3/bwa)
				${BWA} mem -t 6 ${ERV2ref} ${FILENAME_file}.r1.fq ${FILENAME_file}.r2.fq > ${FILENAME_file}.aligned.sam
				echo "BWA has run for ${FILENAME_file}"
				
				#Running samtools to convert sam file to bam file
				${samtools} view -@ 6 -b ${FILENAME_file}.aligned.sam > ${FILENAME_file}.aligned.bam
				
				#Running samtools to sort the files
				${samtools} sort -@ 6 -o ${FILENAME_file}.aligned.sorted.bam ${FILENAME_file}.aligned.bam

				#Running samtools to index the files
				${samtools} index ${FILENAME_file}.aligned.sorted.bam
				echo "Samtools (sam-bam conversion, sort and index) has run for ${FILENAME_file}"
			fi
			fi
	done
echo "END OF ANALYSIS OF BTLTR"

echo "END OF THE SCRIPT"



#######HOW TO RUN THE SCRIPT:
### bash Script_ERVs_Step8_1IlluminaAnalysis_v2.sh > outputfile.txt 2> errorsfile.txt
