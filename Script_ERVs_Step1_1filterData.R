#!/usr/bin/env

#######################################################################
#            PROJECT: Endogenous retroviruses in Bos Taurus           #
#       Subtask: Filtering the data of three cattle populations       #
#             Research Summer School in Statistical Omics             #
#                     author: Alicia L. Bruzos                        #
#                      date: 12th August 2016                         #
#######################################################################

#Set working directory. In my case:
setwd("~/__2016_myProject")

##Step 1: Read files available as a data.frame (separated by tabs)
BBB <- read.table("dataBBB.txt", sep="\t", header=TRUE)
DAM <- read.table("dataDAMONA.txt", sep="\t", header=TRUE)
LIC <- read.table("dataLIC.txt", sep="\t", header=TRUE)

##Step 2: Subset the data which is unique
BBB$DAM <- !(BBB$Breakpoints %in% DAM$Breakpoints)


make the vector column 







##Step 2: Create two new columns in each dataframe comparing the breakpoints with all the breakpoints of other dataframe
#The probability of having two equal breakpoints positions in different chromosomes is almost none
CompareByTwo <- function(df1, df2, df3){
  for i in 1:length(df1$Breakpoints){
      df1$comp_df2 <-  
          if((df1$Breakpoints[i] == df2$Breakpoints[]) || (df1$Breakpoints[] == df2$Breakpoints[])){
                #write 1
          }else{
                #write 0
          }
      df1$comp_df3 <-
          if(df1$Breakpoints[] == df2$Breakpoints[]){
                #write 1
          }else{
                #write 0
          }
  }
}

PopAll <- function(df1, df2, df3){
  if ((df1$comp_df2 == 1) && (df1$comp_df3 == 1)){
      #write in df1$pop <- "All"
  }else if ((df1$comp_df2 == 1) && (df1$comp_df3 == 0)){
      #write in df1$pop <- "df1" && "df2"
  }else if ((df1$comp_df2 == 0) && (df1$comp_df3 == 1)){
      #write in df1$pop <- "df1" && "df3"
  }else{
      #write in df1$pop <- "df1" 
  }}


 






# #Get the ERVs common in the three populations 
# #BBB_DAM_LIC <- BBB[match(BBB$CHROM.START.END, DAM$CHROM.START.END, LIC$CHROM.START.END),]
# 
# #Get the ERVs common only in two populations   (restarle el primero)
# c2_BBB_DAM <- BBB[complete.cases(match(BBB$CHROM.START.END, DAM$CHROM.START.END)),]
# c2_BBB_LIC <- BBB[complete.cases(match(BBB$CHROM.START.END, LIC$CHROM.START.END)),]
# c2_DAM_LIC <- DAM[complete.cases(match(DAM$CHROM.START.END, LIC$CHROM.START.END)),]
# 
# BBB_uniq <- BBB[complete.cases(match(BBB$CHROM.START.END, DAM$CHROM.START.END, LIC$CHROM.START.END), nomatch),]
# DAM_uniq <-
# LIC_uniq <- 
# 
