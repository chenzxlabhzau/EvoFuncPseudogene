library(dplyr)
library(magrittr)

# Expressed dynamic pseudogene
{
  rm(list = ls())
  for (Species in c("Human","Mouse")) {
  wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
  rpkm= read.csv(file.path(wd,Species,paste0(Species,".pseudogene.fpkm.txt")), header = T,sep = "\t")
  coldata <- read.csv(file = file.path("~/Pseudo/Result",Species,"Savedata","coldata.csv"),
                      row.names = 1, sep = ",")
  rpkm = rpkm[,coldata$name]
  if (Species == "Human") {
  fpkm.tissue.stage.max = apply(rpkm, 1, 
                                   function(x)tapply(x, paste0(coldata$condition,".",coldata$stage2),
                                                     max)) %>% t() %>% as.data.frame()
  number.expressed.dynamic.tissues = data.frame()
  for (Tissue in c("Brain","Cerebellum","Heart","Kidney","Liver","Ovary","Testis")) {
    fpkm.tissue.stage.max.choose = dplyr::select(fpkm.tissue.stage.max,starts_with(Tissue))
    stage = colnames(fpkm.tissue.stage.max.choose) %>% 
      lapply(., function(x)strsplit(x,split = ".",fixed = T)[[1]][2]) %>% unlist() %>% as.data.frame()
    stage$time = lapply(stage$.,function(x){substr(x,1,str_length(x)-3)}) %>% 
      unlist() %>% as.numeric()
    stage$unit = lapply(stage$.,function(x){substr(x,str_length(x)-2,str_length(x))}) %>% unlist()
    ddg = read.csv(file = file.path("~/Pseudo/Result",Species,"Savedata/DDG",
                                    paste0(Tissue,"_ddg.csv")))
    if (all(stage$.== paste0(stage$time,stage$unit))) {
      for (i in 1:nrow(stage)) {
        stage[i,"number"] = sum(fpkm.tissue.stage.max.choose[ddg$X,i] >1)
        if (stage[i,3]== "mpb") {
          stage[i,2] = stage[i,2] *30
        }
      }
      if (Tissue != "Ovary") {
        stage$unit %<>% gsub("mpb","dpb",.) %>% factor(.,levels = c("wpc","dpb","ypb"))
      }else{
        stage$unit %<>% factor(.,levels = c("wpc"))
      }
      stage$tissue = Tissue
      number.expressed.dynamic.tissues = rbind(number.expressed.dynamic.tissues,stage)
    }
  }
  number.expressed.dynamic.tissues %<>% arrange(.,unit,time)
  number.expressed.dynamic.tissues$. %<>% factor(.,levels = unique(.))
  p3= number.expressed.dynamic.tissues%>% filter(.,tissue!="Testis")%>% ggplot(.,aes(.,number,color=tissue))+geom_point()+
    geom_smooth(se = TRUE, aes(group = "" ), level = 0.95) + theme_classic()+
    facet_wrap(.~tissue,nrow = 1)+xlab("Developmental stage")+ylab("Expressed dynamic\npseudogene")+
    theme(axis.title.x=element_text(size = 12),legend.position="none",
          axis.text.x = element_text(size=12,angle =45, hjust = 1, vjust = 1),
          axis.title.y=element_text(size = 14),axis.text.y = element_text(size=12),
          strip.background  = element_blank(),
          strip.text = element_text(size=12))+
    scale_x_discrete(breaks=unique(number.expressed.dynamic.tissues$.),
                     labels=c(rep("",6),"10wpc",rep("",7),"0dpb",rep("",19),
                              "14ypb",rep("",12),"58ypb"))+
    scale_color_manual(values=c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397"))
  ggsave(p3,filename = file.path("~/Pseudo/Result",Species,"Picture","Expressed.dynamicpseudogene.noTestis.pdf"),
         device = "pdf",width = 9.5,height = 3)
  p4= number.expressed.dynamic.tissues%>% filter(.,tissue=="Testis")%>% ggplot(.,aes(.,number,color=tissue))+geom_point()+
    geom_smooth(se = TRUE, aes(group = "" ), level = 0.95) + theme_classic()+
    facet_wrap(.~tissue,nrow = 1)+xlab("Developmental stage")+
    theme(axis.title.x=element_text(size = 12),legend.position="none",
          axis.text.x = element_text(size=12,angle =45, hjust = 1, vjust = 1),
          axis.title.y=element_blank(),axis.text.y = element_text(size=12),
          strip.background  = element_blank(),
          strip.text = element_text(size=12))+
    scale_x_discrete(breaks=unique(number.expressed.dynamic.tissues[number.expressed.dynamic.tissues$tissue=="Testis",1]),
                     labels=c(rep("",6),"10wpc",rep("",6),"6mpb",rep("",3),
                              "14ypb",rep("",6),"55ypb"))+
    scale_color_manual(values=c("#ff6700"))
  ggsave(p4,filename = file.path("~/Pseudo/Result",Species,"Picture","Expressed.dynamicpseudogene.Testis.pdf"),
         device = "pdf",width = 5,height = 3)
  
  number.expressed.dynamic.tissues%>% ggplot(.,aes(.,number))+geom_point()+
    geom_smooth(se = TRUE, aes(group = "" ), level = 0.95) + theme_classic()+
    facet_wrap(.~tissue)+ylab("Count")+
    theme(
      strip.background  = element_blank(),
      strip.text = element_text(size=12))
  
  }else if (Species == "Mouse") {
    fpkm.tissue.stage.max = apply(rpkm, 1, 
                                     function(x)tapply(x, paste0(coldata$condition,"_",coldata$stage),
                                                       max)) %>% t() %>% as.data.frame()
    number.expressed.dynamic.tissues = data.frame()
    for (Tissue in c("Brain","Cerebellum","Heart","Kidney","Liver","Ovary","Testis")) {
      fpkm.tissue.stage.max.choose = dplyr::select(fpkm.tissue.stage.max,starts_with(Tissue))
      stage = colnames(fpkm.tissue.stage.max.choose) %>% 
        lapply(., function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist() %>% as.data.frame()
      ddg = read.csv(file = file.path("~/Pseudo/Result",Species,"Savedata/DDG",
                                      paste0(Tissue,"_ddg.csv")))
      for (i in 1:nrow(stage)) {
        stage[i,"number"] = sum(fpkm.tissue.stage.max.choose[ddg$X,i] >1)
        }
      stage$tissue = Tissue
      number.expressed.dynamic.tissues = rbind(number.expressed.dynamic.tissues,stage)
    }
    number.expressed.dynamic.tissues$. %<>% factor(.,levels = unique(coldata$stage))
    p1= number.expressed.dynamic.tissues%>% filter(.,tissue!="Testis")%>% ggplot(.,aes(.,number,color=tissue))+geom_point()+
      geom_smooth(se = TRUE, aes(group = "" ), level = 0.95) + theme_classic()+
      facet_wrap(.~tissue,nrow = 1)+xlab("Developmental stage")+ylab("Expressed dynamic\npseudogene")+
      theme(axis.title.x=element_text(size = 12),legend.position="none",
            axis.text.x = element_text(size=12,angle =45, hjust = 1, vjust = 1),
            axis.title.y=element_text(size = 14),axis.text.y = element_text(size=12),
            strip.background  = element_blank(),
            strip.text = element_text(size=12))+
      scale_x_discrete(breaks=unique(coldata$stage),
                       labels=c("","E11.5","","","","E15.5","","","","0dpb","","","","9wpb"))+
      scale_color_manual(values=c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397","#ff6700"))
    ggsave(p1,filename = file.path("~/Pseudo/Result",Species,"Picture","Expressed.dynamicpseudogene.noTestis.pdf"),
           device = "pdf",width = 9.5,height = 3)
    p2= number.expressed.dynamic.tissues%>% filter(.,tissue=="Testis")%>% ggplot(.,aes(.,number,color=tissue))+geom_point()+
      geom_smooth(se = TRUE, aes(group = "" ), level = 0.95) + theme_classic()+
      facet_wrap(.~tissue,nrow = 1)+xlab("Developmental stage")+ylab("Expressed dynamic\npseudogene")+
      theme(axis.title.x=element_text(size = 12),legend.position="none",
            axis.text.x = element_text(size=12,angle =45, hjust = 1, vjust = 1),
            axis.title.y=element_text(size = 14),axis.text.y = element_text(size=12),
            strip.background  = element_blank(),
            strip.text = element_text(size=12))+
      scale_x_discrete(breaks=unique(coldata$stage),
                       labels=c("","E11.5","","","","E15.5","","","","0dpb","","","","9wpb"))+
      scale_color_manual(values=c("#ff6700"))
    ggsave(p2,filename = file.path("~/Pseudo/Result",Species,"Picture","Expressed.dynamicpseudogene.Testis.pdf"),
           device = "pdf",width = 5,height = 3)
  }
  }
}
