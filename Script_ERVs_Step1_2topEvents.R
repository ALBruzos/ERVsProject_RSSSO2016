#!/usr/bin/env

#######################################################################
#            PROJECT: Endogenous retroviruses in Bos Taurus           #
#           Subtask: Analysis of the Top Events Proportions...        #
#             Research Summer School in Statistical Omics             #
#                     author: Alicia L. Bruzos                        #
#                      date: 12th August 2016                         #
#######################################################################

#Set working directory. In my case:
setwd("~/__2016_myProject")

#Read files available as a data.frame (separated by tabs)
ERVs <- read.table("TopEventsData.txt", sep="\t", header=TRUE)

#Compute the Proportions of each populations and the combined one
ERVs$pBBB <- (ERVs$BBB/sum(ERVs$BBB))*100
ERVs$pDAM <- (ERVs$Damona/sum(ERVs$Damona))*100
ERVs$pLIC <- (ERVs$LIC/sum(ERVs$LIC))*100
ERVs$pCombined <- (ERVs$Combined/sum(ERVs$Combined))*100

#Look for the sorted topEvents for each dataset and the combined one
ERVs_comb <- ERVs[order(ERVs$pCombined, decreasing=1),]
ERVs_BBB <- ERVs[order(ERVs$pBBB, decreasing=1),]
ERVs_DAM <- ERVs[order(ERVs$pDAM, decreasing=1),]
ERVs_LIC <- ERVs[order(ERVs$pLIC, decreasing=1),]

#Construct a table with the top 2 values
ERVs_decision <- data.frame(ERVs_comb[1], ERVs_BBB[1], ERVs_DAM[1], ERVs_LIC[1])

#Change the names of the columns
names(ERVs_decision)[1]<-paste("Combined")
names(ERVs_decision)[2]<-paste("BBB")
names(ERVs_decision)[3]<-paste("DAM")
names(ERVs_decision)[4]<-paste("LIC")

#Exporting the results dataframe into a tab delimited text file
write.table(ERVs_decision, "~/__2016_myProject/TopEventsRanking_Unique.txt", row.names = FALSE, sep="\t", eol="\n")

#Exporting the hole tables of each family...
write.table(ERVs_comb, "~/__2016_myProject/ERVs_comb.txt", row.names = FALSE, sep="\t", eol="\n")
write.table(ERVs_BBB, "~/__2016_myProject/ERVs_BBB.txt", row.names = FALSE, sep="\t", eol="\n")
write.table(ERVs_LIC, "~/__2016_myProject/ERVs_LIC.txt", row.names = FALSE, sep="\t", eol="\n")
write.table(ERVs_DAM, "~/__2016_myProject/ERVs_DAM.txt", row.names = FALSE, sep="\t", eol="\n")



#Still to do the orientation of the ERVs...