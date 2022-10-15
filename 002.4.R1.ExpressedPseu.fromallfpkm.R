rm(list = ls());gc();rm(list = ls())
wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
Num = "002.4.R1."


# A rang of FPKM values ---------------------------------------------------

num = c() 
freq = c()
df = c()
for (S in c("Human","Mouse")) {
  rpkm = fread(file.path(wd,S,paste0(S,".allgene.fpkm.txt"))) %>% as.data.frame()
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  rpkm$type <- b[match(rpkm$V1,b$V4),5]
  rpkm %<>% dplyr::filter(.,grepl("pseudogene",type,ignore.case = T))
  row.names(rpkm) = rpkm$V1
  rpkm = rpkm[,-c(1,ncol(rpkm))]
  rpkm$max = apply(rpkm,1,max)
  c = matrix(0,nrow = 20,ncol = 3) %>% as.data.frame()
  c$V1 = S
  c$V2 = seq(0.1,2,0.1)
  for (i in 1:20) {
    c[i,3] = sum(rpkm$max >= c[i,2]) / nrow(rpkm)
  }
  df = rbind(df,c)
  
  ##################### 
  express = data.frame(gene=row.names(rpkm),max=rpkm$max)
  #Expression level ~ pseudo type
  express$type = b[match(express$gene,b$V4),5]
  express$type2 = express$type %>% gsub("transcribed_","",.) %>% gsub("translated_","",.) %>%
    gsub(" ","",.) %>%
    gsub("unitary_pseudogene","Unitary",.) %>% 
    gsub("polymorphic_pseudogene","Polymorphic",.) %>% 
    gsub("processed_pseudogene","Processed",.)
  express[!express$type2 %in% c("Unitary","Polymorphic","Processed"),"type2"]="Unprocessed"
  print(table(express$type2))
  tmp = matrix(0,nrow = 100 * 4,ncol = 4) %>% as.data.frame()
  tmp$V1 = S
  tmp$V2 = rep(seq(0.1,10,0.1),each=4)
  tmp$V3 = rep(c("Unitary","Polymorphic","Processed","Unprocessed"),100)
  for (i in 1:nrow(tmp)) {
    cutoff = tmp[i,2]; type = tmp[i,3]
    tmp[i,4] = sum(express[express$type2 == type, "max"] >= cutoff) / sum(express$type2 == type)
  }
  freq = rbind(freq,tmp)
}
p1 = ggplot(df,aes(V2,V3))+geom_point(aes(color=V1))+#geom_line()+
  theme_classic()+facet_grid(.~V1)+
  xlab("Expression cutoff (FPKM)")+ylab("(%) Expressed pseudogenes")+
  theme(axis.title.x=element_text(size=12),legend.position='none',
        axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=12),axis.text.y = element_text(size=10),
        strip.background  = element_blank(),
        strip.text = element_text(size=12))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,0.8),
                     breaks = seq(0,0.8,0.2),labels = seq(0,0.8,0.2)*100)
#ggsave(filename = file.path("~/Pseudo/Result","Human","Picture",paste0(Num,"HM.FreqExpressed.pdf")),
#       device = "pdf",height = 3.5,width = 4)
p1
ggplot(df,aes(V2,V3))+geom_point(aes(color=V1))+#geom_line()+
  theme_classic()+#facet_grid(.~V1)+
  xlab("Expression cutoff (FPKM)")+ylab("(%) Expressed pseudogenes")+
  theme(axis.title.x=element_text(size=12),legend.position=c(0.2,0.2),
        legend.title = element_blank(),
        axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=12),axis.text.y = element_text(size=10),
        strip.background  = element_blank(),
        strip.text = element_text(size=12))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,0.8),
                     breaks = seq(0,0.8,0.2),labels = seq(0,0.8,0.2)*100)
ggsave(filename = file.path("~/Pseudo/Result","Human","Picture",paste0(Num,"HM.FreqExpressed.pdf")),
       device = "pdf",height = 3.5,width = 3)

p2 = ggplot(freq,aes(V2,V4))+geom_point(aes(color=V3))+geom_line(aes(group=V3,color=V3))+
  theme_classic()+facet_grid(.~V1)+
  xlab("Expression cutoff (FPKM)")+ylab("(%) Expressed pseudogenes")+
  theme(axis.title.x=element_text(size=12),#legend.position='none',
        axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=12),axis.text.y = element_text(size=10),
        strip.background  = element_blank(),
        strip.text = element_text(size=12))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,1),
                     breaks = seq(0,1,0.1),labels = seq(0,1,0.1)*100)+
  scale_color_manual(values=c("#89c1e7","#5da495","#f7b6ae","#dfcbe2"))+
  guides(color=guide_legend(title="Type"))
ggsave(p2,filename = file.path("~/Pseudo/Result","Human","Picture",paste0(Num,"HM.FreqExpressed.type.pdf")),
       device = "pdf",height = 3.5,width = 6)


# FPKM = 1 as cutoff ---------------------------------------------------
rm(list = ls());gc();rm(list = ls())
wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
Num = "002.4.R1."

num = c() 
freq = c()
for (S in c("Human","Mouse")) {
  rpkm = fread(file.path(wd,S,paste0(S,".allgene.fpkm.txt"))) %>% as.data.frame()
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  rpkm$type <- b[match(rpkm$V1,b$V4),5]
  rpkm %<>% dplyr::filter(.,grepl("pseudogene",type,ignore.case = T))
  row.names(rpkm) = rpkm$V1
  rpkm = rpkm[,-c(1,ncol(rpkm))]
  name = colnames(rpkm) %>% gsub("KidneyTestis","Kidney",.) %>% 
    gsub("Forebrain","Brain",.) %>% 
    gsub("Hindbrain","Cerebellum",.) %>% as.data.frame() 
  #Tissue Stage Sex
  name$TSS = lapply(name$.,function(x)strsplit(as.character(x),"_rep",fixed=T)[[1]][1]) %>% unlist()
  rpkm2 = t(apply(rpkm, 1,  function(x) tapply(x, name$TSS, max))) %>% as.data.frame()
  express = data.frame(gene=row.names(rpkm2),num=apply(rpkm2,1,function(x)sum(x >= 1)))
  write.table(express,file = file.path("~/Pseudo/Result",S,"Savedata","gene.expressnum.1allfpkm.csv"),sep = "\t",quote = F,row.names = F,col.names = T)
}


# Compare expression level between ISO-seq detected and other pseudogenes --------

## cd ~/Pseudo/Data/Ref/Mouse
## grep -i pseudo Mouse.gene.bed | sed 's/chr//g' | bedtools intersect -a - -b PacBio.transcript.bed -wa | sort | uniq > Mouse.PacBiodetected.pseudogene.bed

rm(list = ls());gc();rm(list = ls())
wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
Num = "002.4.R1."

S = "Mouse"
rpkm = fread(file.path(wd,S,paste0(S,".allgene.fpkm.txt"))) %>% as.data.frame()
b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
              header = FALSE,sep = "\t",stringsAsFactors = F)
rpkm$type <- b[match(rpkm$V1,b$V4),5]
rpkm %<>% dplyr::filter(.,grepl("pseudogene",type,ignore.case = T))
row.names(rpkm) = rpkm$V1
rpkm[1:3,340:ncol(rpkm)]
rpkm = rpkm[,-ncol(rpkm)]

rpkm$max = apply(rpkm[,-1],1,max)
pacbio = fread(paste0("~/Pseudo/Data/Seqdata/PacBio/mouse_pseudogene.IsoseqDetected.bed"),header = FALSE) %>% as.data.frame()
rpkm$type = "No"
rpkm[rpkm$V1 %in% pacbio$V3,"type"] = "Yes"
rpkm = rpkm[,c("type","max")]
rpkm[pacbio[grepl("unitary",pacbio$V2),3],]

my_comparisons <- list(c("Yes","No"))
rpkm$type %<>% factor(.,levels = c("Yes","No"))
ggplot(rpkm,aes(type,log10(max+1)))+geom_boxplot(aes(fill=type),notch = TRUE,outlier.colour = "white")+
  theme_classic()+
  xlab("Detected by PacBio")+ylab("Expression level (log10[FPKM])")+
  theme(axis.title.x=element_text(size=14),#legend.position='none',
        axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=14),axis.text.y = element_text(size=12),
        strip.background  = element_blank(),legend.position = "none",
        strip.text = element_text(size=14))+
  stat_compare_means(comparisons = my_comparisons,label.y = 2.75,tip.length = 0)+
  coord_cartesian(ylim = c(0,3))+
  scale_fill_manual(values = c("#B00237","#FAC3D2"))
ggsave(filename = file.path("~/Pseudo/Result",S,"Picture",paste0(Num,"Expressionlevel.IsoseqIllumina.pseudogene.pdf")),
       device = "pdf",width = 2.5,height = 4)

table(rpkm$max >= 1, rpkm$type) 
sum(rpkm[rpkm$type=="Yes",2]>=1) / length(rpkm[rpkm$type=="Yes",2])
sum(rpkm[rpkm$type=="No",2]>=1) / length(rpkm[rpkm$type=="No",2])
#ENSMUSG00000085642

tapply(rpkm$max, rpkm$type, mean)
tapply(rpkm$max, rpkm$type, median)
  