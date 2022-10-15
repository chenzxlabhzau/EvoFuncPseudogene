cd ~/Pseudo/Data/Seqdata/ChuanyuLongqi2019SciDataATAC
grep -v "#" ATAC.address.fromBGI.txt | awk '{print "wget -c ftp://ftp.cngb.org/pub/CNSA/data1/"$1"/"$2"/"$3"/"$4"/*.fq.gz"}' > download.sh
