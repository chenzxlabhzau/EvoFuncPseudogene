
#R
##Heatmap
rm(list = ls())
S="Human"
a <- read.csv("~/Pseudo/Data/Seqdata/WangKuster2019MSBRNAseq/counts.tsv",header = T,sep = "\t",stringsAsFactors = F)
a$X %<>% lapply(.,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist() 
### Choose pseudogene
load(file = file.path("~/Pseudo/Result/GTEx/Savedata/GTEx.RData"))
a$type = m[match(a$X,m$Gene),"type2"]
a %<>% na.omit(.)

#exon length
length <- read.csv(file = file.path("~/Pseudo/Data/Ref",S,paste0(S,".nonreExonlength.bed")),
                   header = FALSE,sep = "\t",stringsAsFactors = FALSE)
#nonredundant gene length
gene.length = length %>% group_by(V4) %>% summarise(length= sum(V5)) %>% as.data.frame()
a$length = gene.length[match(a$X,gene.length$V4),2]


tpm.calculate = function(exprset,len){
  readperlength = t(do.call(rbind, lapply(1:ncol(exprset), function(i){
    exprset[,i]/len})))
  totalcounts <- colSums(readperlength)
  tpm = t(apply(readperlength, 1, function(x) 10^6 * x/totalcounts)) %>% as.data.frame()
  colnames(tpm) = colnames(exprset)
  return(tpm)
}
tpm = tpm.calculate(a[,-c(1,ncol(a)-1,ncol(a))],a$length)

### tissue type
tissue.type = colnames(tpm) %>% lapply(.,function(x)strsplit(x,"_rep",fixed=T)[[1]][1]) %>% unlist() 
### merge by tissue type
tpm = t(apply(tpm, 1, function(x) tapply(x, tissue.type, mean))) %>% as.data.frame()
### mean duplicate genes
tpm$Gene = a$X
tpm = apply(tpm[,-ncol(tpm)], 2, function(x) tapply(x, tpm$Gene, mean)) %>% as.data.frame()

tpm <- tpm[apply(tpm, 1, function(x)sum(x) >0),]
tpm$type = m[match(row.names(tpm),m$Gene),"type2"]
tpm= tpm[order(tpm$type),]
row =  tpm$type %>% as.data.frame()
row.names(row) = row.names(tpm)
colnames(row) ="type"
row$type = factor(row$type,levels = rev(c("Pro","Unp","Poly","Uni")))

save(tpm,row, file = file.path("~/Pseudo/Result/WangKuster2019MSBRNAseq/Savedata/tpm.RData"))

tpm = t(scale(t(tpm[,-ncol(tpm)]))) %>% as.data.frame()
tpm[tpm>2]=2
tpm[tpm<=-2]=-2
ann_colors = list(type=c(Pro="#5da495",Unp="#dfcbe2",Uni="#f7b6ae",Poly="#89c1e7"))
p=pheatmap(tpm,annotation_colors = ann_colors,
           annotation_row = row,
           annotation_legend = T,
           border_color = 'black',#设定每个格子边框的颜色，border=F则无边框
           cluster_rows = F, #对行聚类
           cluster_cols = T, #队列聚类
           show_colnames = T, #是否显示列名
           show_rownames = F, #是否显示行名
           main = "Human (From Wang et al)"
)
ggsave(p,filename = file.path("~/Pseudo/Result/WangKuster2019MSBRNAseq/Picture","Heatmap.tissue.pdf"),device = "pdf",width = 9,height = 6)

## Bar + scatter plot
rm(list = ls())
load(file = file.path("~/Pseudo/Result/WangKuster2019MSBRNAseq/Savedata/tpm.RData"))
tpm <- tpm [,-ncol(tpm)]
tpm$median <- apply(tpm,1,median)
c <- c()
for (i in 1:(ncol(tpm)-1)) {
  c[i] = sum(apply(tpm, 1, function(x) x[i] >x[ncol(tpm)]))
}
colnames(tpm)[which.min(c)]

b= data.frame(c,colnames(tpm)[-ncol(tpm)])
colnames(b)=c("number","tissue")
b= b[order(b$number),]
b$tissue = factor(b$tissue,levels = b$tissue)
ggplot(b,aes(tissue,number))+geom_bar(stat = "identity",position = "dodge",fill=rep(c("#E39E3E","#C7533B"),c(nrow(b)-1,1)))+theme_classic()+
  ylab("NO. higher expressed \npseudogenes")+
  theme(axis.title.x=element_blank(),legend.position="bottom",
        axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
  scale_y_continuous(expand = c(0,0),limits = c(0,10000))
ggsave(filename = file.path("~/Pseudo/Result/WangKuster2019MSBRNAseq/Picture","Number.HigherexpressedPseu.tissue.pdf"),device = "pdf",width = 9.5,height = 6)

p1=ggplot(tpm,aes(log10(median+1),log10(testis+1)))+geom_point(alpha = 0.1)+
  coord_cartesian(xlim = c(0,4),ylim = c(0,4))+theme_bw()+
  xlab("Mean abundance (All Tissues)")+ylab("Mean abundance (Testis)")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  geom_abline(slope=1, intercept=0,colour="#6195C8")+
  annotate("text",x=0.55,y=3.8,size=6,colour="#ff6700",
           label=paste(sum(tpm$testis > tpm$median),
                       paste0(round(sum(tpm$testis > tpm$median)*100 / nrow(tpm)),"%)"),sep = "("))
p1
ggsave(p1,filename = file.path("~/Pseudo/Result/WangKuster2019MSBRNAseq/Picture","Scatter.testis.pdf"),device = "pdf",width = 5,height = 5)

p2=ggplot(tpm,aes(log10(median+1),log10(heart+1)))+geom_point(alpha = 0.1)+
  coord_cartesian(xlim = c(0,4),ylim = c(0,4))+theme_bw()+
  xlab("Mean abundance (All Tissues)")+ylab("Mean abundance (Heart)")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  geom_abline(slope=1, intercept=0,colour="#6195C8")+
  annotate("text",x=3.45,y=0.2,size=6,colour="#c60202",
           label=paste(sum(tpm$heart > tpm$median),
                       paste0(round(sum(tpm$heart < tpm$median)*100 / nrow(tpm)),"%)"),sep = "("))
p2
ggsave(p2,filename = file.path("~/Pseudo/Result/WangKuster2019MSBRNAseq/Picture","Scatter.heart.pdf"),device = "pdf",width = 5,height = 5)

p3=ggplot(tpm,aes(log10(median+1),log10(pancreas+1)))+geom_point(alpha = 0.1)+
  coord_cartesian(xlim = c(0,4),ylim = c(0,4))+theme_bw()+
  xlab("Mean abundance (All Tissues)")+ylab("Mean abundance (Pancreas)")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  geom_abline(slope=1, intercept=0,colour="#6195C8")+
  annotate("text",x=3.45,y=0.2,size=6,colour="#862461",
           label=paste(sum(tpm$pancreas < tpm$median),
                       paste0(round(sum(tpm$pancreas < tpm$median)*100 / nrow(tpm)),"%)"),sep = "("))
p3
ggsave(p3,filename = file.path("~/Pseudo/Result/WangKuster2019MSBRNAseq/Picture","Scatter.pancreas.pdf"),device = "pdf",width = 5,height = 5)