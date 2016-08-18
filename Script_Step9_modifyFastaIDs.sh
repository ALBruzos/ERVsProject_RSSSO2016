#!/bin/bash
#######################################################################
#            PROJECT: Endogenous retroviruses in Bos Taurus           #
#                       Subtask: modifying fasta IDs                  #
#             Research Summer School in Statistical Omics             #
#                     author: Alicia L. Bruzos                        #
#                      date: 18th August 2016                         #
#######################################################################


##modify fasta IDs to make them smaller
#example
#>chr7_105606395.discordant.bam.ERV2.aligned.sorted.bam_ERV2_Full:1-6935	
#output desired
#>chr7_105606395

#Definitive command:
awk '{if (substr($0,0,1) == ">"){ split($0,arra,"."); print ">test_"substr(arra[1],2,20); } else print $0}' test.fa > test_output_cd.fa
#ERV2
awk '{if (substr($0,0,1) == ">"){ split($0,arra,"."); print ">ERV2_BBB_"substr(arra[1],2,20); } else print $0}' ConsensusERV2_BBB_alicia.fa > ConsensusERV2_BBB_def_alicia.fa
awk '{if (substr($0,0,1) == ">"){ split($0,arra,"."); print ">ERV2_LIC_"substr(arra[1],2,20); } else print $0}' ConsensusERV2_LIC_gaston.fa > ConsensusERV2_LIC_def_gaston.fa
awk '{if (substr($0,0,1) == ">"){ split($0,arra,"."); print ">ERV2_DAM_"substr(arra[1],2,20); } else print $0}' ConsensusERV2_DAM_alicia.fa > ConsensusERV2_DAM_def_alicia.fa
#BTL
awk '{if (substr($0,0,1) == ">"){ split($0,arra,"."); print ">BTL_BBB_"substr(arra[1],2,20); } else print $0}' ConsensusBTL_BBB_alicia.fa > ConsensusBTL_BBB_def_alicia.fa
awk '{if (substr($0,0,1) == ">"){ split($0,arra,"."); print ">BTL_LIC_"substr(arra[1],2,20); } else print $0}' ConsensusBTL_LIC_.fa > ConsensusBTL_LIC_def.fa
awk '{if (substr($0,0,1) == ">"){ split($0,arra,"."); print ">BTL_DAM_"substr(arra[1],2,20); } else print $0}' ConsensusBTL_DAM_nad.fa > ConsensusBTL_DAM_def.fa

#other option (random numbers for chromosomes IDs)
awk '/^>/{print ">ERV2_" ++i; next}{print}' test.fa > test_short.fa
##ERV2
awk '/^>/{print ">BBB_ERV2_" ++i; next}{print}' ConsensusERV2_BBB_alicia.fa > ConsensusERV2_BBB_short_alicia.fa
awk '/^>/{print ">LIC_ERV2_" ++i; next}{print}' ConsensusERV2_LIC_gaston.fa > ConsensusERV2_LIC_short_gaston.fa
awk '/^>/{print ">DAM_ERV2_" ++i; next}{print}' ConsensusERV2_DAM_alicia.fa > ConsensusERV2_DAM_short_alicia.fa
##BTL
awk '/^>/{print ">LIC_ERV2_" ++i; next}{print}' ConsensusERV2_LIC_gaston.fa > ConsensusERV2_LIC_short.fa











##counting number of sequences in a fasta file:
grep -c "^>" file.fa

##add something to end of all header lines:
sed 's/>.*/&WHATEVERYOUWANT/' file.fa > outfile.fa

##extract IDs:
grep -o -E "^>\w+" file.fasta | tr -d ">"


