#######################################################################
#            PROJECT: Endogenous retroviruses in Bos Taurus           #
# Subtask: Getting two positions lists (ERV2 and BTLR) for individuals# 
#                     subsets (done in Step3)                         #
#             Research Summer School in Statistical Omics             #
#                     author: Alicia L. Bruzos                        #
#                      date: 16th August 2016                         #
#######################################################################

#Set working directory. In my case:
setwd("~/__2016_myProject/Step6_PositionsFQ_map_1bamPerRegion_IGV(soloORfull)/")

##Step 1: Read files created in Step 3.2 (separated by tabs)
indv26718718_df <- read.table("~/__2016_myProject/Step3_subsets2LIC/indv26718718_LIC.txt", sep="\t", header=TRUE)
indv20205474_df <- read.table("~/__2016_myProject/Step3_subsets2LIC/indv20205474_LIC.txt", sep="\t", header=TRUE)

#Step 2: subsetting the dataframe for both individuals for ERV2 and BTLTR
ERV2_26718718_df <- subset(indv26718718_df, indv26718718_df$Top_Event == "ERV2-1-LTR_BT-ERVK")
BTL_26718718_df <- subset(indv26718718_df, indv26718718_df$Top_Event == "BTLTR1-ERVK")

ERV2_20205474_df <- subset(indv20205474_df, indv20205474_df$Top_Event == "ERV2-1-LTR_BT-ERVK")
BTL_20205474_df <- subset(indv20205474_df, indv20205474_df$Top_Event == "BTLTR1-ERVK")

##Step 3: Look for the Chr:start-end for each one
ERV2_26718718_PosList <- data.frame(ERV2_26718718_df$CHROM.START.END)
BTL_26718718_PosList <- data.frame(BTL_26718718_df$CHROM.START.END)

ERV2_20205474_PosList <- data.frame(ERV2_20205474_df$CHROM.START.END)
BTL_20205474_PosList <- data.frame(BTL_20205474_df$CHROM.START.END)

##Step 4: change the capital C found in Chr:start-end
ERV2_26718718_PosList$lower <- tolower(ERV2_26718718_PosList$ERV2_26718718_df.CHROM.START.END)
BTL_26718718_PosList$lower <- tolower(BTL_26718718_PosList$BTL_26718718_df.CHROM.START.END)
ERV2_20205474_PosList$lower <- tolower(ERV2_20205474_PosList$ERV2_20205474_df.CHROM.START.END)
BTL_20205474_PosList$lower <- tolower(BTL_20205474_PosList$BTL_20205474_df.CHROM.START.END)

##Step 5: merging the positions lists for ERV2 and BTL and removing duplicates
ERV2_PosList <- data.frame(unique(c(ERV2_26718718_PosList$lower, ERV2_20205474_PosList$lower)))
BTL_PosList <- data.frame(unique(c(BTL_26718718_PosList$lower, BTL_20205474_PosList$lower)))

#Step 6: Writing the output in text files
write.table(ERV2_PosList, "~/__2016_myProject/Step5_PositionsFQ_map_1bamPerRegion_IGV(soloORfull)/ERV2_PosList.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t", eol="\n")
write.table(BTL_PosList, "~/__2016_myProject/Step5_PositionsFQ_map_1bamPerRegion_IGV(soloORfull)/BTL_PosList.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t", eol="\n")



