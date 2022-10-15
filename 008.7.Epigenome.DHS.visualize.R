###### Dynamic distribution

# Bar plot ----------------------------------------------------------------
rm(list = ls())
a = fread("~/Pseudo/Data/Seqdata/Roadmap/DHS/all.hg38.gene.DNase.macs2.narrowPeak",header = F)
a = dplyr::filter(a,grepl("pseudogene|protein_coding|lncRNA",V8,ignore.case = T))
a$V8 %<>% gsub("protein_coding","Coding",.)
a[a$V8!="Coding" & a$V8!="lncRNA","V8"] = "Non-dynamic\npseudogene"

ddg = fread("~/Pseudo/Result/Human/Savedata/DDG/all.ddg.csv",header = FALSE)
a[a$V7 %in% ddg$V1,"V8"] = "Dynamic\npseudogene"
b = a %>% group_by(V7,V8) %>% filter(row_number() == 1) %>% ungroup() %>% as.data.frame()

S="Human"
anno = read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)

c = table(b$V8) %>% as.data.frame()
c$all = c(sum(anno$V5==" protein_coding"),length(unique(ddg$V1)),sum(anno$V5==" lncRNA"),
          dplyr::filter(anno,grepl("pseudogene",V5,ignore.case = T)) %>% nrow()- length(unique(ddg$V1)))
c$ratio = c$Freq/c$all
c$Var1 %<>% factor(.,
                   levels = c("Coding","Dynamic\npseudogene","Non-dynamic\npseudogene","lncRNA"))
ggplot(c,aes(Var1,ratio))+geom_bar(aes(fill=Var1),stat = "identity",position = "dodge")+
  ylab("Accessible proportion")+theme_bw()+ 
  theme(axis.title.x=element_blank(),legend.position='none',
        axis.text.x = element_text(size=12,angle = 45,hjust = 1,vjust = 1),
        axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
  scale_fill_manual(values = c("#EE1F26","#61439A","#4F69B5","#862461"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,1))+
  geom_segment(aes(x=1, y=c[1,4]+0.01, xend=2, yend=c[1,4]+0.01))+
  annotate("text", x=1.6, y=c[1,4]-0.02, label="***",size=5)+
  geom_segment(aes(x=2, y=c[2,4]+0.02, xend=3, yend=c[2,4]+0.02))+
  annotate("text", x=2.5, y=c[2,4]+0.03, label="***",size=5)+
  geom_segment(aes(x=3, y=c[3,4]+0.01, xend=4, yend=c[3,4]+0.01))+
  annotate("text", x=3.4, y=c[3,4]-0.02, label="***",size=5)
ggsave(filename = file.path("~/Pseudo/Result/Roadmap/DHS/Picture/Proportion.number.pdf"),
       device = "pdf",width = 5, height = 5)

# Profile plot ----------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
Num = "008.7."

a = fread("~/Pseudo/Data/Seqdata/Roadmap/DHS/test5.gz") %>% as.data.frame()
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
  xlab("Distance from TSS")+ylab("DHS")+
  theme(legend.title = element_blank(),
        #legend.position=c(0.8,0.85),legend.background = element_blank(),legend.key.size = unit(0.3,"cm"),
        axis.title.x=element_text(size=12),axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=12),axis.text.y = element_text(size=10))+
  scale_color_manual(values =  c("#d95f0d","#9ecae1","#61439A","#4F69B5","#862461"))+
  scale_x_continuous(limits = c(0,2000),expand = c(0,0),breaks = seq(0,2000,500),labels = c("-10,000","-5,000","0","5,000","10,000"))+
  guides(tittle=FALSE)
ggsave(filename = file.path(paste0("~/Pseudo/Result/Roadmap/DHS/Picture/",Num,"Profile.TSS.type5.pdf")),
       device = "pdf",width = 3.7, height = 4)

tmp %>% dplyr::filter(.,name!="Coding") %>%
  ggplot(.,aes(pos,value))+geom_line(aes(color=name))+theme_classic()+
  xlab("Distance from TSS")+ylab("H3K27ac")+
  theme(legend.title = element_blank(),
        #legend.position=c(0.8,0.85),legend.background = element_blank(),legend.key.size = unit(0.3,"cm"),
        axis.title.x=element_text(size=12),axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=12),axis.text.y = element_text(size=10))+
  scale_color_manual(values = c("#61439A","#4F69B5","#862461"))+
  scale_x_continuous(limits = c(0,2000),expand = c(0,0),breaks = seq(0,2000,500),labels = c("-10,000","-5,000","0","5,000","10,000"))+
  guides(tittle=FALSE)
ggsave(filename = file.path(paste0("~/Pseudo/Result/Roadmap/DHS/Picture/",Num,"Profile.TSS.type3.pdf")),
       device = "pdf",width = 3.7, height = 4)
