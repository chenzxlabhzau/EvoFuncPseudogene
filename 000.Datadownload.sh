cd ~/
mkdir Pseudo; cd Pseudo/
mkdir Data Result Script
cd Data/
mkdir Ref Seqdata; cd Ref
mkdir Human Mouse; cd Human
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gen code.v34.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.2wayconspseudos.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip *gz
grep -v "#" Homo_sapiens.GRCh38.98.gtf | awk -F ";" '{print $1"\t"$(NF-1)}'| awk -F "\t" '$3=="gene" {print $1"\t"$4"\t"$5"\t"$9"\t"$10"\t"$7"\t"$2}'| sed 's/gene_id //g'| sed 's/gene_biotype //g'| sed 's/"//g' > geneHuman.bed

cd ~/Pseudo/Data/Ref/Mouse
wget ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz
gunzip Mus_musculus.GRCm38.98.gtf.gz
grep -v "#"  Mus_musculus.GRCm38.98.gtf | awk -F ";" '{print $1"\t"$(NF-1)}'| awk -F "\t" '$3=="gene" {print $1"\t"$4"\t"$5"\t"$9"\t"$10"\t"$7"\t"$2}'| sed 's/gene_id //g'| sed 's/gene_biotype //g'| sed 's/"//g' > geneMouse.bed

#for chicken
cd ~/Pseudo/Data/Ref/Chicken
wget ftp://ftp.ensembl.org/pub/release-98/gtf/gallus_gallus/Gallus_gallus.GRCg6a.98.chr.gtf.gz
for i in {1..28} {30..33} Z W; do wget ftp.ensembl.org/pub/release-98/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.chromosome.${i}.fa.gz; done
gunzip *
cat Gallus_gallus.GRCg6a.dna.chromosome.*.fa > Chicken.fa
rm Gallus_gallus.GRCg6a.dna.chromosome.*fa