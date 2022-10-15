https://www.encodeproject.org/search/?searchTerm=chromhmm&type=Annotation&assembly=mm10&biosample_ontology.classification=tissue&organism.scientific_name=Mus+musculus
```{bash}
cd ~/Pseudo/Data/Seqdata/ENCODE/ChromHMM
grep bed.gz chromhmm.ENCODE.txt > chromhmm.ENCODE.sh
wget -i chromhmm.ENCODE.sh

cd ~/Pseudo/Result/ENCODE/ChromHMM
sed 's/^/chr&/g'  ~/Pseudo/Data/Ref/Mouse/geneMouse.bed > geneMouse.chr.bed
for i in `ls  ~/Pseudo/Data/Seqdata/ENCODE/ChromHMM/*bed.gz`; do zcat $i |bedtools intersect -a geneMouse.chr.bed -b -  -wa -wb |awk -F "\t" '{print $4"\t"$11}' >> gene.state.txt ;done
uniq gene.state.txt | sort | uniq > gene.state.uniq.txt