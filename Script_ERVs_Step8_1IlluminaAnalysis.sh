#!/bin/bash
#######################################################################
#            PROJECT: Endogenous retroviruses in Bos Taurus           #
#                       Subtask: 8.1 ILLUMINA DATA                    #
#             Research Summer School in Statistical Omics             #
#                     author: Alicia L. Bruzos                        #
#                      date: 16th August 2016                         #
#######################################################################
#Repeat the script 6.2 for the illumina data, so add a conversion from bam files to fastq files.

set -e  #Stop the script if it has an error

##Step 1: get all the information that we will need: Positions lists of ERV2, BTLTR, fastq files, BWA, samtools, reference genomes.
BWA=/common/WORK/SCHOOL2016/Data/ERVs/bwa/bwa
samtools=/common/WORK/SCHOOL2016/Data/ERVs/samtools-1.3.1/samtools

ERV2ref=/common/WORK/SCHOOL2016/Data/ERVs/refs/ERV2_Full.fa
BTLref=/common/WORK/SCHOOL2016/Data/ERVs/refs/BTLTR1_ref_v1.fa

BBBreads=/common/WORK/SCHOOL2016/Data/ERVs/BBB/

##Step 2: create a list of the .bam files containing the illumina reads 
ls ${BBBreads}*.bam > bamfiles_BBB.txt
BBBfiles=/common/WORK/SCHOOL2016/students/alicia00/bamfiles_BBB.txt

#Results bamfiles_BBB.txt:
# /common/WORK/SCHOOL2016/Data/ERVs/BBB/chr10_102477536.discordant.bam
# /common/WORK/SCHOOL2016/Data/ERVs/BBB/chr10_33763278.discordant.bam
# /common/WORK/SCHOOL2016/Data/ERVs/BBB/chr10_34445397.discordant.bam
# /common/WORK/SCHOOL2016/Data/ERVs/BBB/chr10_42546471.discordant.bam

##Step 3: Create my new positions lists (only start position!!) for ERV2 and LTR from BBBdataset
#done manually with excell
ST_ERV2_BBB=/common/WORK/SCHOOL2016/students/alicia00/Step7/ST_ERV2_BBB.txt
ST_BTL_BBB=/common/WORK/SCHOOL2016/students/alicia00/Step7/ST_BTL_BBB.txt

##Step 4: Do a for loop which will convert .bam files into .fastq files for every file in bamfiles_BBB.txt
BBBfiles_path=(`cat $BBBfiles`)
for i in ${BBBfiles_path[@]}
	do
		echo $i
		OUTPUTname=`echo ${i} | awk '{n=split($0,end_arr,"/"); print end_arr[n]}'`
		${samtools} sort -n  ${i} | ${samtools} fastq -1 ${OUTPUTname}.r1.fq -2 ${OUTPUTname}.r2.fq -
		#####example
		##/common/WORK/SCHOOL2016/Data/ERVs/samtools-1.3.1/samtools sort -n  /common/WORK/SCHOOL2016/Data/ERVs/BBB/chr1_29159451.discordant.bam | /common/WORK/SCHOOL2016/Data/ERVs/samtools-1.3.1/samtools fastq -1 chr1-291.r1.fq -2 chr1-291.r2.fq -
	done

##Step 5: create two lists (read1 && read2) of the .fq files containing the illumina reads 
fq_files=/common/WORK/SCHOOL2016/students/alicia00/Step7/
ls ${fq_files}*r1.fq > R1_files_BBB.txt
ls ${fq_files}*r2.fq > R2_files_BBB.txt
R1=/common/WORK/SCHOOL2016/students/alicia00/Step7/R1_files_BBB.txt
R2=/common/WORK/SCHOOL2016/students/alicia00/Step7/R2_files_BBB.txt

##Step 6: Do two for loops to align and sort the files for ERV and BTL...
#Create bash arrays for the position lists of ERV2 and BTLTR
ERV2pl=(`cat $ST_ERV2_BBB`) #getting all the lines of that file into a bash array
BTLpl=(`cat $ST_BTL_BBB`)

#For loop for ERV2
for itemA in ${ERV2pl[@]}  #for each item of that array do the following
	do
		echo $itemA
		FILENAME_R1=`grep ${itemA} ${R1}`
		FILENAME_R2=`grep ${itemA} ${R2}`
		##/common/WORK/SCHOOL2016/students/alicia00/Step7/chr4_120445493.discordant.bam.r2.fq
		echo "My input is ${FILENAME}"
		
		if [ "${FILENAME}" == "" ]
			then
				echo "There are not fq files corresponding with the start position ${itemA}"
			
		else if [ "${FILENAME}" != "" ]
			then
				OUTPUT_R1=`echo ${FILENAME_R1} | awk '{n=split($0,end_arr,"/"); print end_arr[n]}'`
				OUTPUT_R2=`echo ${FILENAME_R2} | awk '{n=split($0,end_arr,"/"); print end_arr[n]}'`
				##chr4_120445493.discordant.bam.r2.fq

				#Running BWAMEM to align in that file and getting the sam output file (changing for illumina: https://github.com/lh3/bwa)
				${BWA} mem -t 6 ${ERV2ref} ${OUTPUT_R1} ${OUTPUT_R2} > ${OUTPUT_R1}_R2.sam
				
				#Running samtools to convert sam file to bam file
				$samtools view -@ 6 -b ${OUTPUT_R1}_R2.sam > ${OUTPUT_R1}_R2_aligned.bam

				#Running samtools to sort the files
				$samtools sort -@ 6 -o ${OUTPUT_R1}_R2_aligned.sorted.bam ${OUTPUT_R1}_R2_aligned.bam

				#Running samtools to index the files
				$samtools index ${OUTPUT_R1}_R2_aligned.sorted.bam
			fi
			fi
	done

#For loop for BTLTR
for itemB in ${BTLpl[@]}  #for each item of that array do the following
	do
		echo $itemA
		FILENAME_R1=`grep ${itemB} ${R1}`
		FILENAME_R2=`grep ${itemB} ${R2}`
		##/common/WORK/SCHOOL2016/students/alicia00/Step7/chr4_120445493.discordant.bam.r2.fq
		echo "My input is ${FILENAME}"
		
		if [ "${FILENAME}" == "" ]
			then
				echo "There are not fq files corresponding with the start position ${itemB}"
			
		else if [ "${FILENAME}" != "" ]
			then
				OUTPUT_R1=`echo ${FILENAME_R1} | awk '{n=split($0,end_arr,"/"); print end_arr[n]}'`
				OUTPUT_R2=`echo ${FILENAME_R2} | awk '{n=split($0,end_arr,"/"); print end_arr[n]}'`
				##chr4_120445493.discordant.bam.r2.fq

				#Running BWAMEM to align in that file and getting the sam output file (changing for illumina: https://github.com/lh3/bwa)
				${BWA} mem -t 6 ${BTLref} ${OUTPUT_R1} ${OUTPUT_R2} > ${OUTPUT_R1}_R2.sam
				
				#Running samtools to convert sam file to bam file
				$samtools view -@ 6 -b ${OUTPUT_R1}_R2.sam > ${OUTPUT_R1}_R2_aligned.bam

				#Running samtools to sort the files
				$samtools sort -@ 6 -o ${OUTPUT_R1}_R2_aligned.sorted.bam ${OUTPUT_R1}_R2_aligned.bam

				#Running samtools to index the files
				$samtools index ${OUTPUT_R1}_R2_aligned.sorted.bam
			fi
			fi
	done


#Print something to know that the running arrived here
echo "End of the Script 8.1 (illumina alignmets) for BBB dataset"
