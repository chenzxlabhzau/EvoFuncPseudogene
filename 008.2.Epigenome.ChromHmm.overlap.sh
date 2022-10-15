cd ~/Pseudo/Result/Roadmap/Chromhmm
for i in `ls ~/Pseudo/Data/Seqdata/Roadmap/Chromhmm/*bed`; do awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$3-$2"\t"$4}' $i ; done >>all.mnemonics.bed
sed 's/chr//g' all.mnemonics.bed |bedtools intersect -a - -b ~/Pseudo/Data/Ref/Human/geneHuman.bed -wo > all.mnemonics.gene.bed

cd ~/Pseudo/Result/Human/Savedata/DDG
for i in `ls *_ddg.csv`;do awk -F "," '{print $1}' $i| grep ENS; done >> all.ddg.csv
cd ~/Pseudo/Result/Mouse/Savedata/DDG
for i in `ls *_ddg.csv`;do awk -F "," '{print $1}' $i| grep ENS; done >> all.ddg.csv