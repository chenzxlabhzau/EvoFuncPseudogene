cd ~/Pseudo/Data/Seqdata/RPFdb/Human
for i in `cat Human.riboseq.metadata.txt`;do wget http://sysbio.gzzoc.com/rpfdb/data/RPKM/Human_${i}_RPKM.metatable;done

cd ~/Pseudo/Data/Seqdata/RPFdb/Mouse
for i in `cat Mouse.riboseq.metadata.txt`; do wget http://sysbio.gzzoc.com/rpfdb/data/RPKM/Mouse_${i}_RPKM.metatable;done