wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all_hg38lift.mnemonics.bedFiles.tgz

nohup wget https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation.tar.gz &
  
for i in `grep DNase.macs2.narrowPeak.gz DNase.metadata.txt| awk -F "\t" '{print $2}'`; do  echo "nohup wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$i &"; done > down.DNase.sh
bash down.DNase.sh

#18 state
cd ~/Pseudo/Data/Seqdata/Roadmap/Chromhmm/18state
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/all_hg38lift.mnemonics.bedFiles.tgz
tar -zxvf all_hg38lift.mnemonics.bedFiles.tgz  
zcat *bed.gz >> all_18_core_K27ac_hg38lift_mnemonics.bed
bedtools intersect -a ~/Pseudo/Data/Seqdata/Roadmap/Chromhmm/18state/all_18_core_K27ac_hg38lift_mnemonics.bed -b ~/MamDC/Data/Ref/Human/Human.promoter.coding.bed -wa -wb > ~/MamDC/Result/Roadmap/Chromhmm/Savedata/18state/all_18_core_K27ac_hg38lift_mnemonics.promoter.bed

#50 state
cd ~/Pseudo/Data/Seqdata/Roadmap/Chromhmm/50state
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/class1Models_50states/E003/E003_50_segments.bed.gz
nohup wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/class1Models_50states/E004/E004_50_segments.bed.gz &
nohup wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/class1Models_50states/E005/E005_50_segments.bed.gz &
nohup wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/class1Models_50states/E006/E006_50_segments.bed.gz &
nohup wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/class1Models_50states/E007/E007_50_segments.bed.gz &
nohup wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/class1Models_50states/E008/E008_50_segments.bed.gz &
nohup wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/class1Models_50states/E017/E017_50_segments.bed.gz &
for i in `ls *gz`; do zcat $i | /home/qian/source/liftOver stdin ~/Pseudo/Data/Ref/Human/hg19ToHg38.over.chain ${i%bed.gz}hg38.bed unmap; done
cat *hg38.bed >> all.50_segments.hg38.bed
bedtools intersect -a ~/Pseudo/Data/Seqdata/Roadmap/Chromhmm/50state/all.50_segments.hg38.bed -b ~/MamDC/Data/Ref/Human/Human.promoter.coding.bed -wa -wb > ~/MamDC/Result/Roadmap/Chromhmm/Savedata/50state/all.50_segments.hg38.promoter.bed
