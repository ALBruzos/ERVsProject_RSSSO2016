#! /bin/bash
# This script is used to schedule PBS array jobs. Settings for PBS:
### Set the shell PBS should use
#PBS -S /bin/bash
### Queue name
#PBS -q school
### Number of nodes/processors/RAM
#PBS -l select=1:ncpus=6:mem=10gb

#######################################################################
#            PROJECT: Endogenous retroviruses in Bos Taurus           #
#  Subtask: ERV2 alignment and analysis for dataset PB_LIC_p1_2016.fq #
#             Research Summer School in Statistical Omics             #
#                     author: Alicia L. Bruzos                        #
#                      date: 15th August 2016                         #
#######################################################################

#Defining the paths to the programs and datasets:
BWA=/common/WORK/SCHOOL2016/Data/ERVs/bwa/bwa
samtools=/common/WORK/SCHOOL2016/Data/ERVs/samtools-1.3.1/samtools

ERVfullRef=/common/WORK/SCHOOL2016/Data/ERVs/refs/ERV2_Full.fa
p1data=/common/WORK/SCHOOL2016/Data/ERVs/pacbio/PB_LIC_p1_2016.fq

#Step 1: BWAMEM. Running the alignment
$BWA mem -t 6 -x pacbio $ERVfullRef $p1data > PB_p1_alicia.sam

#Step 2: Converting the sam file to bam
$samtools view -@ 6 -b PB_p1_alicia.sam > PB_p1_alicia.bam

#Step 3: Sorting with samtools
$samtools sort -@ 6 -T /common/tmp/PB_p1_alicia.sorted -o PB_p1_alicia.sorted.bam PB_p1_alicia.bam

#Printing something if the script ends well.
echo "Done."