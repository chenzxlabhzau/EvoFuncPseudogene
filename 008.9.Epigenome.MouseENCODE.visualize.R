rm(list = ls());gc();rm(list = ls())
Num = "008.9."
gene.state = fread("~/Pseudo/Result/ENCODE/ChromHMM/gene.state.uniq.txt",header = F)
gene.number = gene.state %>% dplyr::group_by(V1) %>% summarise(number = n()) %>% ungroup() %>% as.data.frame()

gene.anno = fread("~/Pseudo/Result/ENCODE/ChromHMM/geneMouse.chr.bed",header = F)  %>% as.data.frame()

gene.number$type = gene.anno[match(gene.number$V1,gene.anno$V4),5]
gene.number %<>% dplyr::filter(.,grepl("pseudogene|protein_coding|lncRNA",type,ignore.case = T))
gene.number$type %<>% gsub("protein_coding","Coding",.) %>% gsub("lncRNA","lncRNA",.)
gene.number[gene.number$type!="Coding" & gene.number$type!="lncRNA","type"] = "Non-dynamic\npseudogene"
S="Mouse"
ddg = fread(file.path("~/Pseudo/Result",S,"Savedata/DDG/all.ddg.csv"),header = FALSE)
gene.number[gene.number$V1 %in% ddg$V1,"type"] = "Dynamic\npseudogene"

gene.number %<>% dplyr::filter(.,type!="lncRNA")
ggplot(gene.number,aes(number))+geom_bar(aes(y = ..prop..,group=type,fill=type))+
  xlab("Number of ChromHMM states")+ylab("Proportion")+theme_bw()+ 
  theme(axis.title.x=element_text(size=14),legend.position='none',
        axis.text.x = element_text(size=12),
        axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
  scale_fill_manual(values = c("#EE1F26","#61439A","#4F69B5","#862461"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,.85))
ggsave(filename = file.path("~/Pseudo/Result/ENCODE/ChromHMM/Picture",paste0(Num,S,".proportion.number.pdf")),
       device = "pdf",width = 4.5,height = 4)
