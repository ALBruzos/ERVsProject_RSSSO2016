#!/bin/bash
#######################################################################
#            PROJECT: Endogenous retroviruses in Bos Taurus           #
#      Putting all the carriers in one cell instead of several        #
#             Research Summer School in Statistical Omics             #
#                     author: Alicia L. Bruzos                        #
#                      date: 11th August 2016                         #
#######################################################################

#imagine that your starting column of things is the 8th
gawk '{OFS=".\t"; for(i=1; i<=7; i++){print $i}; OFS=";"; for(i=8;i<=$NF;i++){print$i};}' input.txt > output.txt