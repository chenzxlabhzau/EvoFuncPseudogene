rm(list = ls());gc();rm(list = ls())
wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
Num = "002.4."

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
  express = data.frame(gene=row.names(rpkm2),num=apply(rpkm2,1,function(x)sum(x >= 0.3)))
  write.table(express,file = file.path("~/Pseudo/Result",S,"Savedata","gene.expressnum.0.3allfpkm.csv"),sep = "\t",quote = F,row.names = F,col.names = T)
  express.number <- function(expreset,expresscutoff,tissuestagesexcutoff){
    expreset$num = apply(expreset,1,function(x)sum(x >expresscutoff))
    sum(expreset$num >=tissuestagesexcutoff)/nrow(expreset)
  }
  c = matrix(0,nrow = 3,ncol = 3) %>% as.data.frame()
  c$V1 =S
  c$V2 =c("(1,0.1)","(1,0.3)","(3,0.3)") 
  c$V3 =c(express.number(rpkm2,0.1,1),express.number(rpkm2,0.3,1),express.number(rpkm2,0.3,3))
  num = rbind(num,c)
  
  ##################### 
  #Expression level ~ pseudo type
  express$type = b[match(express$gene,b$V4),5]
  express$type2 = express$type %>% gsub("transcribed_","",.) %>% gsub("translated_","",.) %>%
    gsub(" ","",.) %>%
    gsub("unitary_pseudogene","Unitary",.) %>% 
    gsub("polymorphic_pseudogene","Polymorphic",.) %>% 
    gsub("processed_pseudogene","Processed",.)
  express[!express$type2 %in% c("Unitary","Polymorphic","Processed"),"type2"]="Unprocessed"
  #express[grep("unitary",express$type2,ignore.case = T),"type2"]="Unitary"
  table(express$type,express$type2)
  tapply(express$num>0, express$type2, sum)/table(express$type2)
  d = as.data.frame(tapply(express$num>0, express$type2, sum)/table(express$type2))
  d$species = S
  freq = rbind(freq,d)
  
  rpkm$type = express[match(row.names(rpkm),express$gene),"type2"]
  mbg = t(apply(rpkm[,-ncol(rpkm)], 2, function(x)tapply(x, rpkm$type, max))) %>% as.data.frame()
  mbg$tissue = row.names(mbg)
  mbg.longer= pivot_longer(mbg,cols=1:4)
  ggplot(mbg.longer,aes(name,log10(value),fill=name))+geom_boxplot(outlier.color = "white",notch = T)+
    theme_classic()+ylab("Expression level (log10[rpkm+1])")+
    theme(axis.title.x=element_blank(),legend.position='none',
          axis.text.x = element_text(size=14,angle = 45, hjust = 1, vjust = 1),
          axis.title.y=element_text(size=16),axis.text.y = element_text(size=18))+
    scale_fill_manual(values=c("#89c1e7","#5da495","#f7b6ae","#dfcbe2"))
  ggsave(filename = file.path("~/Pseudo/Result",S,"Picture",paste0(S,".Expressionlevel.pseudotype.pdf")),device = "pdf",
         width = 5,height = 4.5)
  ####################
}
p1=ggplot(num,aes(V2,V3))+geom_bar(aes(fill=V1),stat = "identity",position ="dodge")+
  theme_classic()+facet_grid(.~V1)+ylab("Ratio of expressed pseudogenes (%)")+
  theme(axis.title.x=element_blank(),legend.position='none',
        axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=12),axis.text.y = element_text(size=10),
        strip.background  = element_blank(),
        strip.text = element_text(size=12))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,0.8),
                     breaks = seq(0,0.8,0.2),labels = seq(0,0.8,0.2)*100)
#+scale_fill_manual(values=c("#e8ba31","#56a4d8"))
ggsave(p1,filename = file.path("~/Pseudo/Result","Human","Picture",paste0(Num,"HM.FreqExpressed.0.3.pdf")),
       device = "pdf",height = 3.5,width = 3.5)

# Human (1,0.1) 0.7254002
# Human (1,0.3) 0.4897665
# Human (3,0.3) 0.3765416
# Mouse (1,0.1) 0.7341763
# Mouse (1,0.3) 0.5274083
# Mouse (3,0.3) 0.4030113

freq$Freq = freq$Freq * 100
p2 = ggplot(freq,aes(Var1,Freq))+geom_bar(aes(fill=Var1),stat = "identity",position ="dodge")+
  theme_classic()+facet_grid(.~species)+ylab("Ratio of expressed pseudogenes (%)")+
  theme(axis.title.x=element_blank(),legend.position='none',
        axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=12),axis.text.y = element_text(size=10),
        strip.background  = element_blank(),
        strip.text = element_text(size=12))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,100))+
  scale_fill_manual(values = c("#89c1e7","#5da495","#f7b6ae","#dfcbe2"))
ggsave(p2,filename = file.path("~/Pseudo/Result","Human","Picture",paste0(Num,"HM.FreqExpressed.0.3.type.pdf")),
       device = "pdf",height = 3.5,width = 5)

# Polymorphic 45.23810   Human
#   Processed 47.46319   Human
#     Unitary 69.29825   Human
# Unprocessed 51.68174   Human
# Polymorphic 38.63636   Mouse
#   Processed 57.83647   Mouse
#     Unitary 79.72973   Mouse
# Unprocessed 36.21974   Mouse
