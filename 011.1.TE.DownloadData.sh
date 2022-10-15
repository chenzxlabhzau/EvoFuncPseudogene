cd ~/Pseudo/Data/Ref/Human/TE
#grep -E "LINE|SINE|LTR|DNA|Simple_repeat|Low_complexity" Human.hg38UCSC |  awk -F "\t" '{print $6"\t"$7"\t"$8"\t"$12"\t"$13"\t"$10}' > Human.hg38UCSC.repeatmasker.bed
zcat Human.hg38.repeakmasker.UCSC.gz | grep -E "LINE|SINE|LTR|DNA" |  awk -F "\t" '{print $6"\t"$7"\t"$8"\t"$12"\t"$13"\t"$10}' > Human.hg38UCSC.repeatmasker.bed
zcat Mouse.mm10.repeakmasker.UCSC.gz | grep -E "LINE|SINE|LTR|DNA" |  awk -F "\t" '{print $6"\t"$7"\t"$8"\t"$12"\t"$13"\t"$10}' > Mouse.mm10UCSC.repeatmasker.bed

for i in LINE SINE LTR DNA; do grep $i Human.hg38UCSC.repeatmasker.bed > Human.hg38UCSC.${i}.bed; done 
gzip Human.hg38UCSC.*.bed

for i in LINE SINE LTR DNA; do grep $i  Mouse.mm10UCSC.repeatmasker.bed > Mouse.mm10UCSC.${i}.bed; done
gzip Mouse.mm10UCSC*bed