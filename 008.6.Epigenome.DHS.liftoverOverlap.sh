cd ~/Pseudo/Data/Seqdata/Roadmap/DHS
awk -F "\t" '{print $1"\t"$2"\t"$3}' *Peak > all.hg19.DNase.macs2.narrowPeak
/home/qian/source/liftOver all.hg19.DNase.macs2.narrowPeak ~/Pseudo/Data/Ref/Human/hg19ToHg38.over.chain all.hg38.DNase.macs2.narrowPeak   unmap

sed 's/chr//g' all.hg38.DNase.macs2.narrowPeak | bedtools intersect -a - -b ~/Pseudo/Data/Ref/Human/geneHuman.bed -wo > all.hg38.gene.DNase.macs2.narrowPeak

#For profile (window method)

sed 's/^/chr/g' ~/Pseudo/Data/Ref/Human/Genomic.Window.bed | bedtools coverage -a all.hg38.DNase.macs2.narrowPeak -b - > all.hg38.DNase.macs2.window.narrowPeak
awk '{print $1"\t"$2"\t"$3"\t"$7}' all.hg38.DNase.macs2.window.narrowPeak | sort -k1,1 -k2,2n > all.hg38.DNase.macs2.window.sorted.bedGraph
/home/qian/source/bedGraphToBigWig all.hg38.DNase.macs2.window.sorted.bedGraph ~/Pseudo/Data/Ref/Human/Human.chrlength all.hg38.DNase.macs2.window.sorted.bw
computeMatrix reference-point -S all.hg38.DNase.macs2.window.sorted.bw -R ~/Pseudo/Data/Ref/Human/Human.TSS.Coding.bed ~/Pseudo/Data/Ref/Human/Human.TSS.lncRNA.bed ~/Pseudo/Data/Ref/Human/Human.TSS.DynamicPseu.bed ~/Pseudo/Data/Ref/Human/Human.TSS.NodynamicPseu.bed  ~/Pseudo/Data/Ref/Human/Human.TSS.Random.bed -b 10000 -a 10000 -o test5.gz
plotProfile -m test5.gz -o test5.pdf --plotFileFormat pdf
