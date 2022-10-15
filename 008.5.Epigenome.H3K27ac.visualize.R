#####################
#Dynamic distribution

# Bar plot ----------------------------------------------------------------
rm(list = ls())
a = fread("~/Pseudo/Data/Seqdata/Roadmap/H3K27ac/all.H3K27ac.hg38.gene.narrowPeak",
          header = F)
a = dplyr::filter(a,grepl("pseudogene|protein_coding|lncRNA",V9,ignore.case = T))
a$V9 %<>% gsub("protein_coding","Coding",.)
a[a$V9!="Coding" & a$V9!="lncRNA","V9"] = "Non-dynamic\npseudogene"

ddg = fread("~/Pseudo/Result/Human/Savedata/DDG/all.ddg.csv",header = FALSE)
a[a$V8 %in% ddg$V1,"V9"] = "Dynamic\npseudogene"

b= a$V8 %>% unique() %>% as.data.frame()
b$V9 = a[match(b$.,a$V8),9]

S="Human"
anno = read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)

c = table(b$V9) %>% as.data.frame()
c$all = c(sum(anno$V5==" protein_coding"),length(unique(ddg$V1)),sum(anno$V5==" lncRNA"),
          dplyr::filter(anno,grepl("pseudogene",V5,ignore.case = T)) %>% nrow()- length(unique(ddg$V1)))
c$ratio = c$Freq/c$all
c$Var1 %<>% factor(.,
                   levels = c("Coding","Dynamic\npseudogene","Non-dynamic\npseudogene","lncRNA"))
ggplot(c,aes(Var1,ratio))+geom_bar(aes(fill=Var1),stat = "identity",position = "dodge")+
  ylab("Proportion coverred by H3K27ac")+theme_bw()+ 
  theme(axis.title.x=element_blank(),legend.position='none',
        axis.text.x = element_text(size=12,angle = 45,hjust = 1,vjust = 1),
        axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
  scale_fill_manual(values = c("#EE1F26","#61439A","#4F69B5","#862461"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,1))+
  geom_segment(aes(x=1, y=c[1,4]+0.01, xend=2, yend=c[1,4]+0.01))+
  annotate("text", x=1.5, y=c[1,4]+0.015, label="***",size=5)+
  geom_segment(aes(x=2, y=c[2,4]+0.02, xend=3, yend=c[2,4]+0.02))+
  annotate("text", x=2.5, y=c[2,4]+0.03, label="***",size=5)+
  geom_segment(aes(x=3, y=c[3,4]+0.02, xend=4, yend=c[3,4]+0.02))+
  annotate("text", x=3.5, y=c[3,4]+0.03, label="***",size=5)
ggsave(filename = file.path("~/Pseudo/Result/Roadmap/H3K27ac/Picture/Proportion.number.pdf"),
       device = "pdf",width = 5, height = 5)

# profile plot ----------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
Num = "008.5."
if (FALSE) {
  dir = file.path("~/Pseudo/Data/Seqdata/Roadmap/H3K27ac")
  Files = grep("E003-H3K27ac.narrowPeak.gz$",list.files(dir),value=TRUE) 
  #Files = Files[-5]
  filePath <- sapply(Files,function(x){paste(dir,x,sep='/')}) %>% as.list()
  
  S = "Human"
  promoter = fread(file.path("~/Pseudo/Data/Ref",S,paste0(S,".promoter.bed"))) %>% as.data.frame()
  colnames(promoter) = c("chr","start","end","name","strand")
  
  anno = read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                  header = FALSE,sep = "\t",stringsAsFactors = F)
  promoter$type = anno[match(promoter$name,anno$V4),5] %>% gsub(" ","",.)
  promoter %<>% dplyr::filter(.,grepl("pseudogene|protein_coding|lncRNA",type,ignore.case = T))
  promoter$type %<>% gsub("protein_coding","Coding",.)
  promoter[promoter$type!="Coding" & promoter$type!="lncRNA","type"] = "Non-dynamic\npseudogene"
  
  ddg = fread("~/Pseudo/Result/Human/Savedata/DDG/all.ddg.csv",header = FALSE)
  promoter[promoter$name %in% ddg$V1,"type"] = "Dynamic\npseudogene"
  
  random = fread("~/Pseudo/Data/Ref/Human/Human.intergenic.Promoter.frompseudo.seed100.bed") %>% as.data.frame()
  colnames(random) = colnames(promoter)[1:5]
  random$type = "Shuffled"
  
  library(GenomicRanges)
  library(ChIPseeker)
  promoter %<>% rbind(.,random) %>% na.omit %>% GRanges()
  promoter = promoter[promoter@ranges@width==3001,]
  
  pro1 = promoter[promoter$type == "Coding",]
  pro2 = promoter[promoter$type == "Dynamic\npseudogene",]
  pro3 = promoter[promoter$type == "Non-dynamic\npseudogene",]
  pro4 = promoter[promoter$type == "Shuffled",]
  
  tagMatrixList1 <- lapply(filePath, getTagMatrix, windows=pro1)
  tagMatrixList2 <- lapply(filePath, getTagMatrix, windows=pro2)
  tagMatrixList3 <- lapply(filePath, getTagMatrix, windows=pro3)
  tagMatrixList4 <- lapply(filePath, getTagMatrix, windows=pro4)
  
  tagMatrixList = c(tagMatrixList1,tagMatrixList2,tagMatrixList3,tagMatrixList4)
  
  plotAvgProf(tagMatrixList1, xlim=c(-2000, 1000))+theme_classic()+
    plotAvgProf(tagMatrixList2, xlim=c(-2000, 1000))+theme_classic()+
    plotAvgProf(tagMatrixList3, xlim=c(-2000, 1000))+theme_classic()+
    plotAvgProf(tagMatrixList4, xlim=c(-2000, 1000))+theme_classic()
  
  plotAvgProf(tagMatrixList, xlim=c(-2000, 1000))+theme_classic()
  ylab(paste(type[n],"density",sep = " "))+
    theme(legend.title = element_blank(),axis.title.y=element_text(size=13))
  
}

a = fread("~/Pseudo/Data/Seqdata/Roadmap/H3K27ac/test5.gz") %>% as.data.frame()
colnames(a)[1] = "V1"
a = a[,-c(1:3,5,6)]

S = "Human"
cod = fread(file = file.path("~/Pseudo/Data/Ref",S,paste0(S,".TSS.Coding.bed")))
cod$id = paste0(cod$V1,":",cod$V2,"-",cod$V3)

lnc = fread(file = file.path("~/Pseudo/Data/Ref",S,paste0(S,".TSS.lncRNA.bed")))
lnc$id = paste0(lnc$V1,":",lnc$V2,"-",lnc$V3)

dyna = fread(file.path("~/Pseudo/Data/Ref",S,paste0(S,".TSS.DynamicPseu.bed")))
dyna$id = paste0(dyna$V1,":",dyna$V2,"-",dyna$V3)

nody = fread(file.path("~/Pseudo/Data/Ref",S,paste0(S,".TSS.NodynamicPseu.bed")))
nody$id = paste0(nody$V1,":",nody$V2,"-",nody$V3)

a$type = "Shuffled"
a[a$V4 %in% cod$id,"type"] = "Coding"
a[a$V4 %in% lnc$id,"type"] = "lncRNA"
a[a$V4 %in% dyna$id,"type"] = "Dynamic"
a[a$V4 %in% nody$id,"type"] = "Non-dynamic"
a[is.na(a)]=0
tmp = apply(a[,-c(1,ncol(a))], 2, function(x)tapply(x, a$type, function(a)mean(a,trim=0.05))) %>% t %>% as.data.frame()
tmp$pos = 1:nrow(tmp)
tmp %<>% pivot_longer(.,cols=1:c(ncol(.)-1)) %>% as.data.frame()
tmp$name %<>% factor(.,levels = c("Coding","lncRNA","Dynamic","Non-dynamic","Shuffled"))

ggplot(tmp,aes(pos,value))+geom_line(aes(color=name))+theme_classic()+
  xlab("Distance from TSS")+ylab("H3K27ac")+
  theme(legend.title = element_blank(),
        #legend.position=c(0.8,0.85),legend.background = element_blank(),legend.key.size = unit(0.3,"cm"),
        axis.title.x=element_text(size=12),axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=12),axis.text.y = element_text(size=10))+
  scale_color_manual(values =  c("#d95f0d","#9ecae1","#61439A","#4F69B5","#862461"))+
  scale_x_continuous(limits = c(0,2000),expand = c(0,0),breaks = seq(0,2000,500),labels = c("-10,000","-5,000","0","5,000","10,000"))+
  guides(tittle=FALSE)
ggsave(filename = file.path(paste0("~/Pseudo/Result/Roadmap/H3K27ac/Picture/",Num,"Profile.TSS.type5.pdf")),
       device = "pdf",width = 3.7, height = 4)

tmp %>% dplyr::filter(.,name!="Coding") %>%
  ggplot(.,aes(pos,value))+geom_line(aes(color=name))+theme_classic()+
  xlab("Distance from TSS")+ylab("DNase I hypersensitivity")+
  theme(legend.title = element_blank(),
        #legend.position=c(0.8,0.85),legend.background = element_blank(),legend.key.size = unit(0.3,"cm"),
        axis.title.x=element_text(size=12),axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=12),axis.text.y = element_text(size=10))+
  scale_color_manual(values = c("#61439A","#4F69B5","#862461"))+
  scale_x_continuous(limits = c(0,2000),expand = c(0,0),breaks = seq(0,2000,500),labels = c("-10,000","-5,000","0","5,000","10,000"))+
  guides(tittle=FALSE)
ggsave(filename = file.path(paste0("~/Pseudo/Result/Roadmap/DHS/Picture/",Num,"Profile.TSS.type3.pdf")),
       device = "pdf",width = 3.7, height = 4)

