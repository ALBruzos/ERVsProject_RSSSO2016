#! /bin/bash
# This script is used to schedule PBS array jobs.

# Settings for PBS:
### Set the shell PBS should use
#PBS -S /bin/bash

### Queue name
#PBS -q school

### Number of nodes/processors/RAM
#PBS -l select=1:ncpus=6:mem=10gb

#######################################################################
#            PROJECT: Endogenous retroviruses in Bos Taurus           #
#                       Subtask: BWA analysis                         #
#             Research Summer School in Statistical Omics             #
#                     author: Alicia L. Bruzos                        #
#                      date: 15th August 2016                         #
#######################################################################
# TASK: PBS script to take the example data and mapped using BWA to the full ERV reference.
# https://github.com/lh3/bwa
# PBS requirements: 6cpu cores and 10RAM
# Reference: ERV2_full
# Data: example
# Software:
#   *BWAMEM: align (https://github.com/lh3/bwa)
# 		bwa mem -x pacbio ref.fa reads.fq > aln.sam
#		more info: /common/WORK/SCHOOL2016/Data/ERVs/bwa/bwa mem -help
# *SAMTOOLS: convert sam>bam  && sort  (http://www.htslib.org/doc/samtools.html)
# 		
# This script runs the BWA and samtools analysis in several folders as a PBS script.

#Defining the paths to the programs:
BWA=/common/WORK/SCHOOL2016/Data/ERVs/bwa/bwa
samtools=/common/WORK/SCHOOL2016/Data/ERVs/samtools-1.3.1/samtools

#######################################
####STEP 1: BWAMEM                    #
####Running the alignment             #
#######################################
##Example 1 (full path in the script)
$BWA mem -t 6 -x pacbio /common/WORK/SCHOOL2016/Data/ERVs/refs/ERV2_Full.fa /common/WORK/SCHOOL2016/Data/ERVs/pacbio/example1.fq > example1_alicia.sam 
##Example 2 (path in the command line, accessed by $1)
#$BWA mem -t 6 -x pacbio /common/WORK/SCHOOL2016/Data/ERVs/refs/ERV2_Full.fa $1 > example1_alicia.sam

#######################################
####STEP 2: SAMTOOLS                  #
####Converting the sam file to bam    #
#######################################
$samtools view -@ 6 -b example1_alicia.sam > example1_alicia.bam

#######################################
####STEP 3: SAMTOOLS                  #
####Sorting                           #
#######################################
$samtools sort -@ 6 -T /common/tmp/example1_alicia.sorted -o example1_alicia.sorted.bam example1_alicia.bam

echo "Done."
echo name of script is $0
echo first argument is $1

##HOW TO RUN IT?
#chmod 777 Script_ERVs_Step4.sh
#qsub -F /common/WORK/SCHOOL2016/Data/ERVs/pacbio/example1.fq Script_ERVs_Step4.sh
#qstat
#ls