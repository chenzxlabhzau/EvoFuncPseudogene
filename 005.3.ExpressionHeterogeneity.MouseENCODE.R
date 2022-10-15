#bash
#cd ~/Pseudo/Data/Seqdata/ENCODE
#wget https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-36025/E-GEOD-36025.sdrf.txt
#grep -iv CshlLong E-GEOD-36025.sdrf.txt | awk -F "\t" '{print $2"\t"$(NF-22)"\t"$(NF-21)"\t"$20}' | grep SRR > sample.srr.link
#for i in `awk -F "\t" '{print $3}'  sample.srr.link`;do echo "nohup wget $i &" ; done > download.sh
#
#module load snakePipes/1.3.0
#RNA-seq -i /data/processing1/qian/Projects/tmp1/ -o /data/processing1/qian/Projects/tmp2/Snakepipes/ --libraryType 2 -j #30 --DAG --bwBinSize 40 mm10

#R
##Heatmap
rm(list = ls())
Num = "005.3."
S="Mouse"
a = fread("~/Pseudo/Data/Seqdata/ENCODE/counts.tsv") %>% as.data.frame()
a$V1 <- lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
              header = FALSE,sep = "\t",stringsAsFactors = F)
a$type = b[match(a$V1,b$V4),5]
a <- dplyr::filter(a,grepl("pseudogene",type))
m1 <- dplyr::filter(a,grepl("polymorphic_pseudogene",type))
m2 <- dplyr::filter(a,grepl("unprocessed_pseudogene",type))
m3 <- dplyr::filter(a,grepl("unitary_pseudogene",type))
m4 <- dplyr::filter(a,V1 %in% setdiff(a$V1,c(m1$V1,m2$V1,m3$V1)))
m <- rbind(m1,m2,m3,m4)
m$type2 = rep(c("Poly","Unp","Uni","Pro"),c(nrow(m1),nrow(m2),nrow(m3),nrow(m4)))

#exon length
length <- read.csv(file = file.path("~/Pseudo/Data/Ref",S,paste0(S,".nonreExonlength.bed")),
                     header = FALSE,sep = "\t",stringsAsFactors = FALSE)
#nonredundant gene length
gene.length = length %>% group_by(V4) %>% summarise(length= sum(V5)) %>% as.data.frame()
m$length = gene.length[match(m$V1,gene.length$V4),2]

tpm.calculate = function(exprset,len){
  readperlength = t(do.call(rbind, lapply(1:ncol(exprset), function(i){
    exprset[,i]/len})))
  totalcounts <- colSums(readperlength)
  tpm = t(apply(readperlength, 1, function(x) 10^6 * x/totalcounts)) %>% as.data.frame()
  return(tpm)
}
tpm = tpm.calculate(m[,-c(1,ncol(m),ncol(m)-1,ncol(m)-2)],m$length)
colnames(tpm) = colnames(m)[-c(1,ncol(m),ncol(m)-1,ncol(m)-2)]
name = colnames(tpm) %>% as.data.frame()

info = fread("~/Pseudo/Data/Seqdata/ENCODE/sample.srr.link",header = F) %>% as.data.frame()
name$V2 = info[match(name$.,info$V2),1] 
name$V2 %<>% gsub("CNS","CentralNervousSystem",.) %>%
  gsub("Bladder","UrinaryBladder",.) %>%  gsub("SubcFatPad","SubcutaneousAdiposeTissue",.) %>% gsub("GenitalFatPad","GenitalAdiposeTissue",.) %>%
  gsub("LgIntestine","LargeIntestine",.) %>% gsub("Adrenal","AdrenalGland",.) %>%
  gsub("SmIntestine","SmallIntestine",.)
if (all(colnames(tpm) == name$.)) {
  mbg = t(apply(tpm, 1, function(x) tapply(x, name$V2, mean))) %>% as.data.frame()
}

row.names(mbg) = m$V1
mbg$type = m$type2
row =  mbg$type %>% as.data.frame()
row.names(row) = row.names(mbg)
colnames(row) ="type"
row$type = factor(row$type,levels = c("Pro","Unp","Poly","Uni"))
save(mbg,row, file = file.path("~/Pseudo/Result/ENCODE/Savedata/tpm.RData"))

mbg2 = t(scale(t(mbg[,-ncol(mbg)]))) %>% as.data.frame()
table(colnames(mbg2)[apply(mbg2, 1, which.max) %>% as.numeric()])
mbg2[mbg2>2]=2
mbg2[mbg2<=-2]=-2
ann_colors = list(type=c(Pro="#5da495",Unp="#dfcbe2",Uni="#f7b6ae",Poly="#89c1e7"))
library(pheatmap)
p=pheatmap(mbg2,annotation_colors = ann_colors,
         annotation_row = row,
         annotation_legend = T,
         border_color = 'black',#设定每个格子边框的颜色，border=F则无边框
         cluster_rows = F, #对行聚类
         cluster_cols = T, #队列聚类
         show_colnames = T, #是否显示列名
         show_rownames = F, #是否显示行名
         main = "Mouse (ENCODE)"
)
ggsave(p,filename = file.path("/home/qians/Pseudo/Result/ENCODE/Picture",paste0(Num,"Heatmap.tissue.pdf")),device = "pdf",width = 9,height = 6)

## bar plot
rm(list = ls())
Num = "005.3."
load(file = file.path("~/Pseudo/Result/ENCODE/Savedata/tpm.RData"))
mbg <- mbg [,-ncol(mbg)]
mbg$median <- apply(mbg,1,median)
c <- c()
for (i in 1:(ncol(mbg)-1)) {
    c[i] = sum(apply(mbg, 1, function(x) x[i] >x[ncol(mbg)]))
}
colnames(mbg)[which.min(c)]
b= data.frame(c,colnames(mbg)[-ncol(mbg)])
colnames(b)=c("number","tissue")
b= b[order(b$number),]
b$tissue %<>% gsub("GenitalAdiposeTissue","GAT",.) %>% 
  gsub("SubcutaneousAdiposeTissue","SAT",.) %>% gsub("CentralNervousSystem","CNS",.) 
b$tissue = factor(b$tissue,levels = b$tissue)
p1 = ggplot(b,aes(tissue,number))+geom_bar(stat = "identity",position = "dodge",
                                      fill=rep(c("#E39E3E","#ff6700","#3298c8","#349a00","#E39E3E","#3298c8"),c(17,1,4,1,1,1)))+
  theme_classic()+ylab("Highly expressed pseudogenes")+
    theme(axis.title.x=element_blank(),legend.position="bottom",
          axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1),
          axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
    scale_y_continuous(expand = c(0,0),limits = c(0,4000))
ggsave(p1,filename = file.path("~/Pseudo/Result/ENCODE/Picture",paste0(Num,"Number.HigherexpressedPseu.tissue.pdf"))
       ,device = "pdf",width = 7,height = 5)
p1
barplot(1:7,col= c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397","#ff6700"))
ggplot(b,aes(tissue,number))+geom_bar(stat = "identity",position = "dodge",fill=rep(c("#E39E3E","#C7533B","#E39E3E","#C7533B"),c(18,4,2,1)))+theme_classic()+ylab("NO. higher expressed \npseudogenes")+
    theme(axis.title.x=element_blank(),legend.position="bottom",
          axis.text.x = element_text(size=8,angle = 90, hjust = 1, vjust = 1),
          axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
    scale_y_continuous(expand = c(0,0),limits = c(0,4000),position = "right")
