############### Mouse
cd ~/Pseudo/Data/Seqdata/GTRDforTF/Mouse
wget http://gtrd.biouml.org/downloads/19.10/chip-seq/Mus%20musculus_macs2_clusters.interval.gz

for i in {1..50} X Y; do grep -w "chr"$i Mus\ musculus_macs2_clusters.interval | awk -F "\t" '{print $1"aaa"$6"\t"$2"\t"$3"\t"$5}' >> Mouse.macs2.clusters.interval.bed; done
bedtools sort -i Mouse.macs2.clusters.interval.bed | bedtools sort -i - > Mouse.macs2.clusters.interval.merge.bed
bedtools sort -i Mouse.macs2.clusters.interval.tmp | bedtools merge -i - | sed 's/aaa/\t/g'| awk -F "\t" '{print $1"\t"$3"\t"$4"\t"$2}' > Mouse.macs2.clusters.interval.merge.bed

cd /home/qians/Pseudo/Data/Ref/Mouse
awk -F "\t" '$3=="exon" {print $0}' Mus_musculus.GRCm38.98.gtf |awk -F ";" '{print $1}'| sed 's/gene_id //g' | sed 's/"//g'| awk -F "\t" '{print $9"chr"$7"chr"$1"\t"$4"\t"$5}'| bedtools sort -i - | bedtools merge -i - | sed 's/chr/\t/g'| awk -F "\t" '{print $3"\t"$4"\t"$5"\t"$1"\t"$2}' > Mouse.nonreExonlength.2.bed

############### Human
cd ~/Pseudo/Data/Seqdata/GTRDforTF/Human
awk -F "\t" '{print $1"aaa"$5"\t"$2"\t"$3}' Homo\ sapiens_macs2_clusters.interval | sed '1d' > Human.macs2.clusters.interval.bed
sort -k 1,1 -k 2,2n Human.macs2.clusters.interval.bed |bedtools merge -i - | sed 's/aaa/\t/g'| awk -F "\t" '{print $1"\t"$3"\t"$4"\t"$2}' > Human.macs2.clusters.interval.merge.bed
awk -F " " '{print  $1"\t"$2"\t"$3"\t"$4}' /home/qians/Pseudo/Data/Seqdata/GTRDforTF/Human/Human.macs2.clusters.interval.merge.bed > /home/qians/Pseudo/Data/Seqdata/GTRDforTF/Human/Human.macs2.clusters.interval.merge2.bed

awk -F "\t" '{print $6"xxx"$1"\t"$2"\t"$3"\t"$5}' Homo\ sapiens_macs2_clusters.interval > Human.macs2.clusters.interval2.bed
sed -i '1d' Human.macs2.clusters.interval2.bed
sort -k 1,1 -k 2,2n Human.macs2.clusters.interval2.bed |bedtools merge -i - | sed 's/xxx/\t/g'| awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$1}' > Human.macs2.clusters.interval.merge3.bed


cd /home/qians/Pseudo/Data/Ref/Human
awk -F "\t" '$3=="exon" {print $0}' Homo_sapiens.GRCh38.98.gtf |awk -F ";" '{print $1}'| sed 's/gene_id //g' | sed 's/"//g'| awk -F "\t" '{print $9"chr"$7"chr"$1"\t"$4"\t"$5}'| bedtools sort -i - | bedtools merge -i - | sed 's/chr/\t/g'| awk -F "\t" '{print $3"\t"$4"\t"$5"\t"$1"\t"$2}' > Human.nonreExonlength.2.bed
