cd ~/Pseudo/Data/Seqdata/Roadmap/H3K27ac
awk -F "\t" '{print $1"\t"$2"\t"$3"\t""xYx"NR}' *narrowPeak > all.H3K27ac.hg19.narrowPeak
/home/qian/source/liftOver all.H3K27ac.hg19.narrowPeak ~/Pseudo/Data/Ref/Human/hg19ToHg38.over.chain all.H3K27ac.hg38.narrowPeak unmap
sed 's/chr//g' all.H3K27ac.hg38.narrowPeak| bedtools intersect -a - -b ~/Pseudo/Data/Ref/Human/geneHuman.bed -wo > all.H3K27ac.hg38.gene.narrowPeak

#For profile (window method)
cd ~/Pseudo/Data/Ref/Human
bedtools makewindows -g Human.length -w 100 -s 100 >  ~/Pseudo/Data/Ref/Human/Genomic.Window.bed
cd ~/Pseudo/Data/Seqdata/Roadmap/H3K27ac
sed 's/^/chr/g' ~/Pseudo/Data/Ref/Human/Genomic.Window.bed | bedtools coverage -a all.H3K27ac.hg38.narrowPeak -b - > all.H3K27ac.hg38.window.narrowPeak

awk '{print $1"\t"$2"\t"$3"\t"$7}' all.H3K27ac.hg38.window.narrowPeak | sort -k1,1 -k2,2n > all.H3K27ac.hg38.window.sorted.bedGraph
/home/qian/source/bedGraphToBigWig all.H3K27ac.hg38.window.sorted.bedGraph ~/Pseudo/Data/Ref/Human/Human.chrlength all.H3K27ac.hg38.window.sorted.bw

#Random TSS as negative control 
grep -vE "KI|GL" geneHuman.bed | bedtools flank -i - -g Human.length -b 100000 | sed 's/^/chr/g' | bedtools shuffle -i Human.TSS.Coding.bed -g Human.chrlength -excl - > Human.TSS.Random.bed
computeMatrix reference-point -S all.H3K27ac.hg38.window.sorted.bw -R ~/Pseudo/Data/Ref/Human/Human.TSS.Coding.bed ~/Pseudo/Data/Ref/Human/Human.TSS.lncRNA.bed ~/Pseudo/Data/Ref/Human/Human.TSS.DynamicPseu.bed ~/Pseudo/Data/Ref/Human/Human.TSS.NodynamicPseu.bed  ~/Pseudo/Data/Ref/Human/Human.TSS.Random.bed -b 10000 -a 10000 -o test5.gz
plotProfile -m test5.gz -o test5.pdf --plotFileFormat pdf