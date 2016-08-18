#!/bin/bash
#######################################################################
#            PROJECT: Endogenous retroviruses in Bos Taurus           #
#                       Subtask: 6.2                                  #
#             Research Summer School in Statistical Omics             #
#                     author: Alicia L. Bruzos                        #
#                      date: 16th August 2016                         #
#######################################################################
# TASK: 
#		1) Take all the fq files from one folder (each fq file = 1 site/position of LIC)
#		2) Get two lists of positions: one for BLTR and one for ERV2. [Rscript 5.1]
#		3) Map them to the reference elements. Use a for loop. If they do not do that continue with the next one.
#		4) Get a bam file per region
#		5) Then look in IGV if it is a full-length of solo-LTR element
#

#set -e  #Stop the script if it has an error 

##Step 1: get all the information that we will need: Positions lists of ERV2, BTLTR, fastq files, BWA, samtools, reference genomes.
#generation of fastq_PosList.txt file by running: ls /common/WORK/SCHOOL2016/Data/ERVs/PBsubs/*.fastq > fastq_PosList.txt
fastq_PL=/common/WORK/SCHOOL2016/students/alicia00/Step5/fastq_PosList.txt
ERV2_PL=/common/WORK/SCHOOL2016/students/alicia00/Step5/PL_erv2.txt
BTL_PL=/common/WORK/SCHOOL2016/students/alicia00/Step5/PL_btl.txt

BWA=/common/WORK/SCHOOL2016/Data/ERVs/bwa/bwa
samtools=/common/WORK/SCHOOL2016/Data/ERVs/samtools-1.3.1/samtools

ERV2ref=/common/WORK/SCHOOL2016/Data/ERVs/refs/ERV2_Full.fa
BTLref=/common/WORK/SCHOOL2016/Data/ERVs/refs/BTLTR1_ref_v1.fa

##Step 2: loop to compare the ERV2 and BTL positions lists with the fastq position list files

#Create bash arrays for the position lists of ERV2 and BTLTR
ERV2pl=(`cat $ERV2_PL`) #getting all the lines of that file into a bash array
BTLpl=(`cat $BTL_PL`)

#${ERV2pl[@]} #accessing and returning every item of that array

#For loop for ERV2
for itemA in ${ERV2pl[@]}  #for each item of that array do sthg
	do
		echo $itemA
		
		FILENAME=`grep ${itemA} ${fastq_PL}`
		##/common/WORK/SCHOOL2016/Data/ERVs/PBsubs/chr10:23850844-23851844.fastq
		echo "My input is ${FILENAME}"
		if [ "${FILENAME}" == "" ]
			then
				echo "There is not a fastq file corresponding with the position ${itemA}"
			
		else if [ "${FILENAME}" != "" ]
			then
				OUTPUT=`echo ${FILENAME} | awk '{n=split($0,end_arr,"/"); print end_arr[n]}'`
				##chr10:23850844-23851844.fastq

				#Running BWAMEM to align in that file and getting the sam output file
				$BWA mem -t 6 -x pacbio ${ERV2ref} ${FILENAME} > ${OUTPUT}.sam

				#Running samtools to convert sam file to bam file
				$samtools view -@ 6 -b ${OUTPUT}.sam > ${OUTPUT}.bam

				#Running samtools to sort the files
				$samtools sort -@ 6 -o ${OUTPUT}.sorted.bam ${OUTPUT}.bam

				#Running samtools to index the files
				$samtools index ${OUTPUT}.sorted.bam
			fi
			fi
	done

#For loop for BTLTR 
for itemB in ${BTLpl[@]}  #for each item of that array do sthg
	do
		FILENAME=`grep ${itemB} ${fastq_PL}`
		##/common/WORK/SCHOOL2016/Data/ERVs/PBsubs/chr10:23850844-23851844.fastq
		
		if [ "${FILENAME}" == "" ]
			then
				echo "There is not a fastq file corresponding with the position ${itemB}"
			
		else if [ "${FILENAME}" != "" ]
			then
				OUTPUT=`echo ${FILENAME} | awk '{n=split($0,end_arr,"/"); print end_arr[n]}'`
				##chr10:23850844-23851844.fastq
				
				#Running BWAMEM to align in that file and getting the sam output file
				$BWA mem -t 6 -x pacbio ${BTLref} ${FILENAME} > ${OUTPUT}.sam
				
				#Running samtools to convert sam file to bam file 
				$samtools view -@ 6 -b ${OUTPUT}.sam > ${OUTPUT}.bam
				
				#Running samtools to sort the files 
				$samtools sort -@ 6 -o ${OUTPUT}.sorted.bam ${OUTPUT}.bam
				
				#Running samtools to index the files 
				$samtools index ${OUTPUT}.sorted.bam
			fi
			fi
	done

#Print something to know that the running arrived here
echo "End of the Script 6.2"
