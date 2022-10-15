rm(list = ls());gc();rm(list = ls())
Num = "S007."

pseu = fread("/opt/qians/Pseudogene/Data/UCSC.phastcons/Human/Human.pseudogene.exon.unique.phastcons.30way.bed") %>% as.data.frame()
pseu$V1 %<>% lapply(.,function(x)strsplit(x,split = ":",fixed = TRUE)[[1]][1]) %>% unlist()
pseu = tapply(pseu$V6, pseu$V1, max) %>% as.data.frame()

# DDP ~ conserved ----------------------------------------------------------
S = "Human"
ddpseu = fread(paste0("~/Pseudo/Result/",S,"/Savedata/DDG/all.ddg.csv"),header = FALSE) %>% as.data.frame()

pseu$dd = "Nondynamic"
pseu[row.names(pseu) %in% ddpseu$V1,"dd"]= "Dynamic"
testP = wilcox.test(pseu[pseu$dd=="Dynamic",1],pseu[pseu$dd!="Dynamic",1])
testP = testP$p.value
ggplot(pseu,aes(dd,.))+geom_boxplot(aes(fill=dd),outlier.color = "white",notch = TRUE)+
  theme_bw()+
  ylab("Conservation")+guides(fill=FALSE)+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=10),axis.title.y=element_text(size=12),
        axis.text.y = element_text(size=10))+
  #scale_y_continuous(expand = c(0, 0))+
  geom_segment(aes(x=1, y=0.98, xend=2, yend=0.98),color="black",size=0.2)+
  annotate("text", x=1.5, y=1.02, label=paste0("P= ",format(testP,digits = 2)),size=3)+
  scale_fill_manual(values = c("#61439A","#4F69B5"))
ggsave(filename = paste0("/home/qians/Pseudo/Result/Human/Picture/",Num,"Human.Conservation.dynamic.pdf"),
       width = 3, height =3)


# translated ~ conserved ----------------------------------------------------------
S = "Human"
translate = fread(paste0("~/Pseudo/Result/RPFdb/",S,"/Savedata/",S,".translated.validated.csv"),header = FALSE) %>% as.data.frame()

pseu$dd = "Nontranslated"
pseu[row.names(pseu) %in% translate$V3,"dd"]= "Translated"
testP = wilcox.test(pseu[pseu$dd=="Translated",1],pseu[pseu$dd!="Translated",1])
testP = testP$p.value
pseu$dd %<>% factor(.,levels = c("Translated","Nontranslated"))
ggplot(pseu,aes(dd,.))+geom_boxplot(aes(fill=dd),outlier.color = "white",notch = TRUE)+
  theme_bw()+
  ylab("Conservation")+guides(fill=FALSE)+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size=10,angle = 45,hjust = 1,vjust = 1),axis.title.y=element_text(size=12),
        axis.text.y = element_text(size=10))+
  #scale_y_continuous(expand = c(0, 0))+
  geom_segment(aes(x=1, y=0.92, xend=2, yend=0.92),color="black",size=0.2)+
  annotate("text", x=1.5, y=0.97, label=paste0("P= ",format(testP,digits = 2)),size=3)+
  scale_fill_manual(values = c("#49a88f","#aadde0"))
ggsave(filename = paste0("/home/qians/Pseudo/Result/Human/Picture/",Num,"Human.Conservation.translated.pdf"),
       width = 2, height = 3.5)
