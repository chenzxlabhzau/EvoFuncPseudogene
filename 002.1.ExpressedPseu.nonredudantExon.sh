cd ~/Pseudo/Data/Ref/Human
awk -F "\t" '$3=="exon" {print $0}' Homo_sapiens.GRCh38.98.gtf|awk -F ";" '{print $1}'| sed 's/gene_id //g' | sed 's/"//g'| awk -F "\t" '{print $9"chr"$1"\t"$4"\t"$5}' | bedtools sort -i - | bedtools merge -i - | sed 's/chr/\t/g'| awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$1"\t"$4-$3}' > Human.nonreExonlength.bed
awk -F "\t" '{print $4}' Human.nonreExonlength.bed|sort| uniq |wc -l

# cd ~/Pseudo/Data/Ref/Mouse
# awk -F "\t" '$3=="exon" {print $0}' Mus_musculus.GRCm38.98.gtf |awk -F ";" '{print $1}'| sed 's/gene_id //g' | sed 's/"//g'| awk -F "\t" '{print $9"chr"$1"\t"$4"\t"$5}' | bedtools sort -i - | bedtools merge -i - | sed 's/chr/\t/g'| awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$1"\t"$4-$3}' > Mouse.nonreExonlength.bed
cd ~/Pseudo/Result/PacBio/Savedata
awk -F "\t" '$3=="exon" {print $0}' mouse_final.gtf | awk -F ";" '{print $1}'| sed 's/gene_id //g' | sed 's/"//g' | awk -F "\t" '{print $9"chr"$1"\t"$4"\t"$5}' | bedtools sort -i - | bedtools merge -i - | sed 's/chr/\t/g'| awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$1"\t"$4-$3}' > Mouse.nonreExonlength.PacBio.bed