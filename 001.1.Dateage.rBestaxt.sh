###########################
### download annotation files
cd /home/qians/Pseudo/Data/Ref/
wget -P Chicken/ ftp://ftp.ensembl.org/pub/release-101/gtf/gallus_gallus/Gallus_gallus.GRCg6a.101.gtf.gz
wget -P Chimp/ ftp://ftp.ensembl.org/pub/release-101/gtf/pan_troglodytes/Pan_troglodytes.Pan_tro_3.0.101.gtf.gz # Pan_tro 3.0 = panTro5 https://www.ncbi.nlm.nih.gov/assembly/GCF_000001515.7/
wget -P Dog/ ftp://ftp.ensembl.org/pub/release-101/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.CanFam3.1.101.gtf.gz #
wget -P Opossum/ ftp://ftp.ensembl.org/pub/release-101/gtf/monodelphis_domestica/Monodelphis_domestica.ASM229v1.101.gtf.gz #
wget -P Platypus/ ftp://ftp.ensembl.org/pub/release-101/gtf/ornithorhynchus_anatinus/Ornithorhynchus_anatinus.mOrnAna1.p.v1.101.gtf.gz
wget -P Rabbit/ ftp://ftp.ensembl.org/pub/release-101/gtf/oryctolagus_cuniculus/Oryctolagus_cuniculus.OryCun2.0.101.gtf.gz
wget -P Rat/ ftp://ftp.ensembl.org/pub/release-101/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.101.gtf.gz
wget -P Rhesus/ ftp://ftp.ensembl.org/pub/release-101/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.101.gtf.gz
wget -P XTropicalis/ ftp://ftp.ensembl.org/pub/release-101/gtf/xenopus_tropicalis/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.101.gtf.gz
wget -P Zebrafish/ ftp://ftp.ensembl.org/pub/release-101/gtf/danio_rerio/Danio_rerio.GRCz11.101.gtf.gz

### generate gene bed file
#intersect to find the homologous type of gene
## generate gene bed file for bedtools
for species in `ls`; do zcat ~/Pseudo/Data/Ref/$species/*gtf.gz| grep -w gene | awk -F ";" '{print $1"\t"$(NF-1)"\t"$(NF-3)}' | awk -F "\t" '{print "chr"$1"\t"$4"\t"$5"\t"$9"\t"$10"\t"$NF}'| sed 's/gene_id //g'| sed 's/gene_biotype //g'| sed 's/gene_name //g'| sed 's/"//g' > ~/Pseudo/Data/Ref/$species/$species.gene.bed; done
grep -w gene ~/Pseudo/Data/Ref/Mouse/Mus_musculus.GRCm38.98.gtf | awk -F ";" '{print $1"\t"$(NF-1)"\t"$(NF-3)}' | awk -F "\t" '{print "chr"$1"\t"$4"\t"$5"\t"$9"\t"$10"\t"$NF}'| sed 's/gene_id //g'| sed 's/gene_biotype //g'| sed 's/gene_name //g'| sed 's/"//g' > ~/Pseudo/Data/Ref/Mouse/Mouse.gene.bed

###########################
### download reciprocal chain files
#for human
cd /opt/qians/Pseudogene/Data/UCSC.maf/Human
nohup wget -c -P Chicken/ https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/vsHg38/reciprocalBest/galGal6.hg38.rbest.chain.gz &
nohup wget -c -P Chicken/ https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/vsHg38/reciprocalBest/hg38.galGal6.rbest.chain.gz &
nohup wget -c -P Chimp/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro5/reciprocalBest/panTro5.hg38.rbest.chain.gz &
nohup wget -c -P Chimp/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro5/reciprocalBest/hg38.panTro5.rbest.chain.gz &
nohup wget -c -P Dog/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsCanFam3/reciprocalBest/canFam3.hg38.rbest.chain.gz  &
nohup wget -c -P Dog/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsCanFam3/reciprocalBest/hg38.canFam3.rbest.chain.gz &
nohup wget -c -P Mouse/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMm10/reciprocalBest/hg38.mm10.rbest.chain.gz &
nohup wget -c -P Mouse/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMm10/reciprocalBest/mm10.hg38.rbest.chain.gz &
nohup wget -c -P Opossum/ https://hgdownload.soe.ucsc.edu/goldenPath/monDom5/vsHg38/reciprocalBest/hg38.monDom5.rbest.chain.gz &
nohup wget -c -P Opossum/ https://hgdownload.soe.ucsc.edu/goldenPath/monDom5/vsHg38/reciprocalBest/monDom5.hg38.rbest.chain.gz &
nohup wget -c -P Platypus/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsOrnAna1/reciprocalBest/hg38.ornAna1.rbest.chain.gz &
nohup wget -c -P Platypus/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsOrnAna1/reciprocalBest/ornAna1.hg38.rbest.chain.gz &
nohup wget -c -P Rabbit/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsOryCun2/reciprocalBest/hg38.oryCun2.rbest.chain.gz &
nohup wget -c -P Rabbit/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsOryCun2/reciprocalBest/oryCun2.hg38.rbest.chain.gz &
nohup wget -c -P Rat/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsRn6/reciprocalBest/hg38.rn6.rbest.chain.gz &
nohup wget -c -P Rat/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsRn6/reciprocalBest/rn6.hg38.rbest.chain.gz &
nohup wget -c -P Rhesus/ https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsHg38/reciprocalBest/hg38.rheMac10.rbest.chain.gz &
nohup wget -c -P Rhesus/ https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsHg38/reciprocalBest/rheMac10.hg38.rbest.chain.gz &
nohup wget -c -P XTropicalis https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsXenTro9/reciprocalBest/hg38.xenTro9.rbest.chain.gz &
nohup wget -c -P XTropicalis https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsXenTro9/reciprocalBest/xenTro9.hg38.rbest.chain.gz &
nohup wget -c -P Zebrafish https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsDanRer11/reciprocalBest/danRer11.hg38.rbest.chain.gz &
nohup wget -c -P Zebrafish https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsDanRer11/reciprocalBest/hg38.danRer11.rbest.net.gz &

#for mouse
cd /opt/qians/Pseudogene/Data/UCSC.maf/Mouse
nohup wget -c -P Chicken/ https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/vsMm10/reciprocalBest/galGal6.mm10.rbest.chain.gz &
nohup wget -c -P Chicken/ https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/vsMm10/reciprocalBest/mm10.galGal6.rbest.chain.gz &
nohup wget -c -P Chimp/ https://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsPanTro5/reciprocalBest/mm10.panTro5.rbest.chain.gz &
nohup wget -c -P Chimp/ https://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsPanTro5/reciprocalBest/panTro5.mm10.rbest.chain.gz &
nohup wget -c -P Dog/ https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsMm10/reciprocalBest/canFam4.mm10.rbest.chain.gz & #neeconvert to canfam3
nohup wget -c -P Dog/ https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsMm10/reciprocalBest/mm10.canFam4.rbest.chain.gz &
nohup ##Human &
nohup wget -c -P Opossum/ http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsMonDom5/reciprocalBest/mm10.monDom5.rbest.chain.gz &
nohup wget -c -P Opossum/ http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsMonDom5/reciprocalBest/monDom5.mm10.rbest.chain.gz &
nohup wget -c -P platypus/ http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsOrnAna1/reciprocalBest/mm10.ornAna1.rbest.chain.gz &
nohup wget -c -P platypus/ http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsOrnAna1/reciprocalBest/ornAna1.mm10.rbest.chain.gz &
nohup wget -c -P Rabbit/ https://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsOryCun2/reciprocalBest/mm10.oryCun2.rbest.chain.gz &
nohup wget -c -P Rabbit/ https://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsOryCun2/reciprocalBest/oryCun2.mm10.rbest.chain.gz &
nohup #wget -c -P Rat/ https://hgdownload.soe.ucsc.edu/goldenPath/rn6/vsMm10/reciprocalBest/mm10.rn6.rbest.chain.gz &
nohup #wget -c -P Rat/ https://hgdownload.soe.ucsc.edu/goldenPath/rn6/vsMm10/reciprocalBest/rn6.mm10.rbest.chain.gz &
nohup wget -c -P Rhesus/ https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsMm10/reciprocalBest/mm10.rheMac10.rbest.chain.gz &
nohup wget -c -P Rhesus/ https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsMm10/reciprocalBest/rheMac10.mm10.rbest.chain.gz &
nohup wget -c -P XTropicalis https://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsXenTro9/reciprocalBest/mm10.xenTro9.rbest.chain.gz &
nohup wget -c -P XTropicalis https://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsXenTro9/reciprocalBest/xenTro9.mm10.rbest.chain.gz &
nohup wget -c -P Zebrafish https://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsDanRer11/reciprocalBest/danRer11.mm10.rbest.chain.gz &
nohup wget -c -P Zebrafish https://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsDanRer11/reciprocalBest/mm10.danRer11.rbest.chain.gz &

###########################
## liftover 
# Human
cd ~/Pseudo/Data/Ref/Human
grep pseu Homo_sapiens.GRCh38.98.gtf | grep -w exon |awk -F ";" '{print $1"\t"$8"\t"$12}' | sed 's/gene_id //g'|sed 's/gene_biotype //g'| sed 's/exon_id //g'| sed 's/"//g'  | awk -F "\t" '{print $1"\t"$4"\t"$5"\t"$9":"$11":"$5-$4}'  | sed 's/ //g' > Human.pseudogene.exon.bed
#Align one species to Human, and align back
cd /opt/qians/Pseudogene/Data/UCSC.maf/Human
for species in `ls`; do chain1=`ls $species/hg38.*.rbest.chain`; chain2=`ls $species/*.hg38.rbest.chain`; echo "# == 1.align human to $species =="; echo "# == 2.re-align $species to human =="; sed 's/^/chr/g' ~/Pseudo/Data/Ref/Human/Human.pseudogene.exon.bed | sort -k1,1 -k2,2n | uniq| /home/qian/source/liftOver stdin $chain1 stdout /dev/null | /home/qian/source/liftOver stdin $chain2 stdout /dev/null | sed 's/:/\t/g' | awk -F "\t" '{print $0"\t"$3-$2}'| awk '$7*2>=$6 {print $0}' > $species/Human.$species.Human.pseudogene.exon.alignment.bed;done
#align satisfied human pseudogenes to one species to find homologous regions
for species in `ls`; do chain1=`ls $species/hg38.*rbest.chain`; echo "# == re-align human pseudogene above 50% to $species to detect ortholog regions =="; awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4":"$5}' $species/Human.$species.Human.pseudogene.exon.alignment.bed | /home/qian/source/liftOver stdin $chain1 $species/Human.$species.Human.$species.pseudogene.exon.alignment.bed /dev/null ; done

# Mouse
cd ~/Pseudo/Data/Ref/Mouse
grep pseu Mus_musculus.GRCm38.98.gtf | grep -w exon |awk -F ";" '{print $1"\t"$8"\t"$12}' | sed 's/gene_id //g'|sed 's/gene_biotype //g'| sed 's/exon_id //g'| sed 's/"//g' | awk -F "\t" '{print $1"\t"$4"\t"$5"\t"$9":"$11":"$5-$4}' | sed 's/ //g' > Mouse.pseudogene.exon.bed
cd /opt/qians/Pseudogene/Data/UCSC.maf/Mouse/Human
ln -s /opt/qians/Pseudogene/Data/UCSC.maf/Human/Mouse/*rbest.chain ./
cd /opt/qians/Pseudogene/Data/UCSC.maf/Mouse
gunzip */*gz
#Align one species to Mouse, and align back
for species in `ls`; do chain1=`ls $species/mm10.*.rbest.chain`; chain2=`ls $species/*.mm10.rbest.chain`; echo "# == 1.align mouse to $species =="; echo "# == 2.re-align $species to mouse =="; sed 's/^/chr/g' ~/Pseudo/Data/Ref/Mouse/Mouse.pseudogene.exon.bed | sort -k1,1 -k2,2n | uniq| /home/qian/source/liftOver stdin $chain1 stdout /dev/null | /home/qian/source/liftOver stdin $chain2 stdout /dev/null | sed 's/:/\t/g' | awk -F "\t" '{print $0"\t"$3-$2}'| awk '$7*2>=$6 {print $0}' > $species/Mouse.$species.Mouse.pseudogene.exon.alignment.bed;done
#align satisfied mouse pseudogenes to one species to find homologous regions
for species in `ls`; do chain1=`ls $species/mm10.*rbest.chain`; echo "# == re-align mouse pseudogene above 50% to $species to detect ortholog regions =="; awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4":"$5}' $species/Mouse.$species.Mouse.pseudogene.exon.alignment.bed | /home/qian/source/liftOver stdin $chain1 $species/Mouse.$species.Mouse.$species.pseudogene.exon.alignment.bed /dev/null ; done

###########################
## use bed file to overlap homologous regions with ref annotated gene
# Human
cd /opt/qians/Pseudogene/Data/UCSC.maf/Human
for species in `ls `; do echo "overlap homologous regions with $species ref annotated gene" ;bedtools intersect -a $species/Human.$species.Human.$species.pseudogene.exon.alignment.bed -b ~/Pseudo/Data/Ref/$species/$species.gene.bed -wa -wb > $species/$species.overlap.Humanexon.bed; done


###########################
## filter paramter: no overlap with TE, no Y
#Human
cd ~/Pseudo/Data/Ref/Human/TE
zcat Human.hg38.repeakmasker.UCSC.gz | awk -F "\t" '{print $6"\t"$7"\t"$8"\t"$12"\t"$13"\t"$10}' | sed '1d' | sed 's/chr//g' > Human.hg38.repeakmasker.UCSC.bed
bedtools intersect -a ~/Pseudo/Data/Ref/Human/Human.pseudogene.exon.bed -b Human.hg38.repeakmasker.UCSC.bed -wo > Human.pseudogene.exon.repeat.txt

cd ~/Pseudo/Data/Ref/Mouse/TE
zcat Mouse.mm10.repeakmasker.UCSC.gz | awk -F "\t" '{print $6"\t"$7"\t"$8"\t"$12"\t"$13"\t"$10}' | sed '1d' | sed 's/chr//g' > Mouse.mm10.repeakmasker.UCSC.bed
bedtools intersect -a ~/Pseudo/Data/Ref/Mouse/Mouse.pseudogene.exon.bed -b Mouse.mm10.repeakmasker.UCSC.bed -wo > Mouse.pseudogene.exon.repeat.txt

#Mouse
