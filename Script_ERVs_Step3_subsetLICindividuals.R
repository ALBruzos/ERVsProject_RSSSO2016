#!/usr/bin/env

#######################################################################
#            PROJECT: Endogenous retroviruses in Bos Taurus           #
#                  Subtask: Subset 2 LIC individuals                  #
#             Research Summer School in Statistical Omics             #
#                     author: Alicia L. Bruzos                        #
#                      date: 13th August 2016                         #
#######################################################################

#Set working directory. In my case:
setwd("~/__2016_myProject/Step0_datasets/")

#######################################################
### LIC dataset subsetted for individual 26718718 and 
### 20205474    
#######################################################

##Step 1: Read files available as a data.frame (separated by tabs)
BBB <- read.table("dataBBB.txt", sep="\t", header=TRUE)
DAM <- read.table("dataDAMONA.txt", sep="\t", header=TRUE)
LIC <- read.table("dataLIC_TopEventsFixed.txt", sep="\t", header=TRUE)

# Individual 26718718
indv26718718_df <- LIC[grep("26718718", LIC),]
write.table(indv26718718_df, "~/__2016_myProject/Step3_subsets&positionsLists/indv26718718_LIC.txt", row.names = FALSE, sep="\t", eol="\n")

# Individual 20205474
indv20205474_df <- LIC[grep("20205474", LIC),]
write.table(indv20205474_df, "~/__2016_myProject/Step3_subsets&positionsLists/indv20205474_LIC.txt", row.names = FALSE, sep="\t", eol="\n")


### RESULTS:
##Individual    Cases
## 26718718     167
## 20205474     304
