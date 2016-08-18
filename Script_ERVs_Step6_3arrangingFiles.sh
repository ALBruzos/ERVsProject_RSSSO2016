#BASH
#arranging the bam files created into two folders to them categorising them using IGV

mkdir BTL_bams
mkdir ERV2_bams

ls *.bam | grep -f PL_erv2.txt | xargs mv -t ERV2_bams/
ls *.bam | grep -f PL_btl.txt | xargs mv -t BTL_bams/

cp *.bai ERV2_bams/
cp *.bai BTL_bams/

rename : _ *
rename - _ *