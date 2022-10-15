#Mouse number
cd ~/Pseudo/Result/GTRDforTF/Mouse
bedtools intersect -a ~/Pseudo/Data/Ref/Mouse/Mouse.promoter.bed -b ~/Pseudo/Data/Seqdata/GTRDforTF/Mouse/Mouse.macs2.clusters.interval.merge.bed  -c > Mouse.promoter.number.bed
###Promoter region shuffle
cd ~/Pseudo/Data/Ref/Mouse
grep -vE "^JH|^GL" geneMouse.bed | bedtools flank -i - -b 1000 -g Mouse.length | bedtools shuffle -i Mouse.promoter.bed -excl - -g Mouse.length -seed 100 > Mouse.intergenic.Promoter.frompseudo.seed100.bed
cd ~/Pseudo/Result/GTRDforTF/Mouse
bedtools intersect -a ~/Pseudo/Data/Ref/Mouse/Mouse.intergenic.Promoter.frompseudo.seed100.bed -b ~/Pseudo/Data/Seqdata/GTRDforTF/Mouse/Mouse.macs2.clusters.interval.merge.bed -c > Mouse.randompromoter.number.bed

#Mouse type
bedtools intersect -a ~/Pseudo/Data/Ref/Mouse/Mouse.promoter.bed -b ~/Pseudo/Data/Seqdata/GTRDforTF/Mouse/Mouse.macs2.clusters.interval.merge.bed -wa -wb > Mouse.promoter.type1.bed
awk -F "\t" '{print $4"NAME"$NF}' Mouse.promoter.type1.bed| sort | uniq > Mouse.promoter.type2.bed
awk -F "NAME" '{print $1}' Mouse.promoter.type2.bed | sort | uniq -c >Mouse.promoter.type3.bed

#Human number
cd ~/Pseudo/Result/GTRDforTF/Human/
bedtools intersect -a ~/Pseudo/Data/Ref/Human/Human.promoter.bed -b ~/Pseudo/Data/Seqdata/GTRDforTF/Human/Human.macs2.clusters.interval.merge2.bed -c > Human.promoter.number.bed
## For random shuffle
cd ~/Pseudo/Data/Ref/Human
###gene region shuffle
bedtools shuffle -i geneHuman.pseudogene.bed -excl geneHuman.bed -g Human.length -seed 100 > Human.intergenic.Gene.frompseudo.seed100.bed #gene region shuffle
###Promoter region shuffle
grep -vE "^KI|^GL" geneHuman.bed | bedtools flank -i - -b 1000 -g Human.length | bedtools shuffle -i  Human.promoter.bed -excl - -g Human.length -seed 100 > Human.intergenic.Promoter.frompseudo.seed100.bed
cd ~/Pseudo/Result/GTRDforTF/Human/
bedtools intersect -a ~/Pseudo/Data/Ref/Human/Human.intergenic.Promoter.frompseudo.seed100.bed -b ~/Pseudo/Data/Seqdata/GTRDforTF/Human/Human.macs2.clusters.interval.merge2.bed -c > Human.randompromoter.number.bed

cd ~/Pseudo/Result/GTRDforTF/Human/
sed 's/^/chr/g' ~/Pseudo/Data/Ref/Human/Human.intergenic.Gene.frompseudo.seed100.bed | bedtools intersect -a - -b ~/Pseudo/Data/Seqdata/GTRDforTF/Human/Human.macs2.clusters.interval.merge.bed -c > Human.intergenic.Gene.frompseudo.seed100.TFnumber.bed 
#Human type
bedtools intersect -a Human.promoter.bed -b Human.macs2.clusters.interval.merge3.bed -wa -wb > Human.promoter.type1.bed
awk -F "\t" '{print $4"NAME"$NF}' Human.promoter.type1.bed | sort | uniq |awk -F "NAME" '{print $1}' | sort | uniq -c > Human.promoter.type3.bed