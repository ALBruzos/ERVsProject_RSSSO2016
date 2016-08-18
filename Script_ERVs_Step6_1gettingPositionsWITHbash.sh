


alicia00@lobsang:~/Step5> grep 20205474 ../LIC_dataToSubset.txt >> LICsIndv.txt
alicia00@lobsang:~/Step5> grep 26718718 ../LIC_dataToSubset.txt >> LICsIndv.txt
alicia00@lobsang:~/Step5> grep ERV2-1-LTR_BT-ERVK LICsIndv.txt | cut -f 1 | awk '{print tolower($0)}' > PL_erv2.txt
alicia00@lobsang:~/Step5> grep BTLTR1-ERVK LICsIndv.txt | cut -f 1 | awk '{print tolower($0)}' > PL_btl.txt
