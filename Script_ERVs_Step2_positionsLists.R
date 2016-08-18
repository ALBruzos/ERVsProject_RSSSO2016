#!/usr/bin/env

#######################################################################
#            PROJECT: Endogenous retroviruses in Bos Taurus           #
#   Subtask: Subset data to analyse unique ERVs to do the phylogeny   #
#             Research Summer School in Statistical Omics             #
#                     author: Alicia L. Bruzos                        #
#                      date: 13th August 2016                         #
#######################################################################

#Set working directory. In my case:
setwd("~/__2016_myProject/Step0_datasets/")

#######################################################
### SUBTASK: 01
### List of positions of all the populations which have
### ERV2-1-LTR_BT-ERVK  and  BTLTR1-ERVK
#######################################################

##Step 1: Read files available as a data.frame (separated by tabs)
BBB <- read.table("dataBBB.txt", sep="\t", header=TRUE)
DAM <- read.table("dataDAMONA.txt", sep="\t", header=TRUE)
LIC <- read.table("dataLIC_TopEventsFixed.txt", sep="\t", header=TRUE)

#Delete from the LIC$Top_Event the numbers after ":"
#LIC$Top_Event_2 <- grepl(":",)

##Step 2: Look for the CHR:start-end for ERV2-1-LTR_BT-ERVK	for the three populations; 
##then save it in a text file
ERV2_BBB_df <- subset(BBB, BBB$Top_Event == "ERV2-1-LTR_BT-ERVK")
ERV2_BBB_PosList <- data.frame(ERV2_BBB_df$CHROM.START.END)
write.table(ERV2_BBB_PosList, "~/__2016_myProject/ERV2_BBB_PosList.txt", row.names = FALSE, sep="\t", eol="\n")

ERV2_LIC_df <- subset(LIC, LIC$Top_Event == "ERV2-1-LTR_BT-ERVK")
ERV2_LIC_PosList <- data.frame(ERV2_LIC_df$CHROM.START.END)
write.table(ERV2_LIC_PosList, "~/__2016_myProject/ERV2_LIC_PosList.txt", row.names = FALSE, sep="\t", eol="\n")

ERV2_DAM_df <- subset(DAM, DAM$Top_Event == "ERV2-1-LTR_BT-ERVK")
ERV2_DAM_PosList <- data.frame(ERV2_DAM_df$CHROM.START.END)
write.table(ERV2_DAM_PosList, "~/__2016_myProject/ERV2_DAM_PosList.txt", row.names = FALSE, sep="\t", eol="\n")


##Step 3: Look for the CHR:start-end for BTLTR1-ERVK
BTL_BBB_df <- subset(BBB, BBB$Top_Event == "BTLTR1-ERVK")
BTL_BBB_PosList <- data.frame(BTL_BBB_df$CHROM.START.END)
write.table(BTL_BBB_PosList, "~/__2016_myProject/BTL_BBB_PosList.txt", row.names = FALSE, sep="\t", eol="\n")

BTL_LIC_df <- subset(LIC, LIC$Top_Event == "BTLTR1-ERVK")
BTL_LIC_PosList <- data.frame(BTL_LIC_df$CHROM.START.END)
write.table(BTL_LIC_PosList, "~/__2016_myProject/BTL_LIC_PosList.txt", row.names = FALSE, sep="\t", eol="\n")

BTL_DAM_df <- subset(DAM, DAM$Top_Event == "BTLTR1-ERVK")
BTL_DAM_PosList <- data.frame(BTL_DAM_df$CHROM.START.END)
write.table(BTL_DAM_PosList, "~/__2016_myProject/BTL_DAM_PosList.txt", row.names = FALSE, sep="\t", eol="\n")


###RESULTS:
##        BTLTR1-ERVK   ERV2-1-LTR_BT-ERVK
##BBB     233           182
##LIC     408           313
##DAM     305           254


#######################################################
### SUBTASK: 02
### LIC dataset subsetted for individual 26718718 and 
### 20205474    
#######################################################

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


