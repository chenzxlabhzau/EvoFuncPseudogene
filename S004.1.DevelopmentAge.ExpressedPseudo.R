rm(list = ls());gc();rm(list = ls())
Num = "S004.1."
library(dplyr)
library(magrittr)

# Fraction of expressed dynamic pseudogene
{
  rm(list = ls())
  Num = "S004.1."
  for (Species in c("Human","Mouse")) {
    wd = "/home/qians/Pseudo/Data/Seqdata/Illumina"
    rpkm= fread(file.path(wd,Species,paste0(Species,".allgene.fpkm.txt"))) %>% as.data.frame()
    row.names(rpkm) = rpkm$V1
    rpkm = rpkm[,-1]
    rpkm %<>% dplyr::select(contains("rep"))
    coldata <- read.csv(file = file.path("~/Pseudo/Result",Species,"Savedata","coldata.csv"),
                        row.names = 1, sep = ",")
    b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",Species,paste0("gene",Species,".bed")),
                  header = FALSE,sep = "\t",stringsAsFactors = F)
    rpkm$type <- b[match(row.names(rpkm),b$V4),5]
    rpkm = rpkm[grepl("pseu",rpkm$type),coldata$name]
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
        
        age = read.csv(file = file.path("~/Pseudo/Result",Species,"Savedata",
                                        paste0(Species,".gene.age.csv")))
        ddg = read.csv(file = file.path("~/Pseudo/Result",Species,"Savedata/DDG",
                                        paste0(Tissue,"_ddg.csv")))
        if (all(stage$.== paste0(stage$time,stage$unit))) {
          for (i in 1:nrow(stage)) {
            stage[i,"number"] = sum(fpkm.tissue.stage.max.choose[ddg$X,i] > 0.3) / sum(fpkm.tissue.stage.max.choose[,i] > 0.3) * 100
              #mean(age[age[,1] %in% row.names(fpkm.tissue.stage.max.choose[fpkm.tissue.stage.max.choose[,i]>0.3,]),"age"])
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
      p3= number.expressed.dynamic.tissues %>% ggplot(.,aes(.,number,color=tissue))+geom_point()+
        geom_smooth(se = TRUE, aes(group = "" ), level = 0.95,method = "lm") + theme_classic()+
        facet_wrap(.~tissue,nrow = 1,scales = "free_x")+xlab("Developmental stage")+ylab("Fraction of expressed dynamic\npseudogene (%)")+
        theme(axis.title.x=element_text(size = 12),legend.position="none",
              axis.text.x = element_text(size=12,angle =45, hjust = 1, vjust = 1),
              axis.title.y=element_text(size = 14),axis.text.y = element_text(size=12),
              strip.background  = element_blank(),
              strip.text = element_text(size=12))+
        scale_x_discrete(breaks=unique(number.expressed.dynamic.tissues$.),
                         labels=c(rep("",6),"10wpc",rep("",7),"0dpb",rep("",19),
                                  "14ypb",rep("",12),"58ypb"))+
        scale_color_manual(values=c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397","#ff6700"))
      ggsave(p3,filename = file.path("~/Pseudo/Result",Species,"Picture",paste0(Num,"Fraction.ExpressedDynamicPseudogene.pdf")),
             device = "pdf",width = 9.5,height = 4.5)
      
    }else if (Species == "Mouse") {
      fpkm.tissue.stage.max = apply(rpkm, 1, 
                                    function(x)tapply(x, paste0(coldata$condition,"_",coldata$stage),
                                                      max)) %>% t() %>% as.data.frame()
      number.expressed.dynamic.tissues = data.frame()
      for (Tissue in c("Brain","Cerebellum","Heart","Kidney","Liver","Ovary","Testis")) {
        fpkm.tissue.stage.max.choose = dplyr::select(fpkm.tissue.stage.max,starts_with(Tissue))
        stage = colnames(fpkm.tissue.stage.max.choose) %>% 
          lapply(., function(x)strsplit(x,split = "_",fixed = T)[[1]][2]) %>% unlist() %>% as.data.frame()
        age = read.csv(file = file.path("~/Pseudo/Result",Species,"Savedata",
                                        paste0(Species,".gene.age.csv")))
        
        ddg = read.csv(file = file.path("~/Pseudo/Result",Species,"Savedata/DDG",
                                        paste0(Tissue,"_ddg.csv")))
        
        for (i in 1:nrow(stage)) {
          stage[i,"number"] = sum(fpkm.tissue.stage.max.choose[ddg$X,i] > 0.3) / sum(fpkm.tissue.stage.max.choose[,i] > 0.3) * 100
            #mean(age[age[,1] %in% row.names(fpkm.tissue.stage.max.choose[fpkm.tissue.stage.max.choose[,i]>0.3,]),"age"])
          
        }
        stage$tissue = Tissue
        number.expressed.dynamic.tissues = rbind(number.expressed.dynamic.tissues,stage)
      }
      number.expressed.dynamic.tissues$. %<>% factor(.,levels = unique(coldata$stage))
      p1 = number.expressed.dynamic.tissues %>% ggplot(.,aes(.,number,color=tissue))+geom_point()+
        geom_smooth(se = TRUE, aes(group = "" ), level = 0.95,method = "lm") + theme_classic()+
        facet_wrap(.~tissue,nrow = 1)+xlab("Developmental stage")+ylab("Fraction of expressed dynamic\npseudogene (%)")+
        theme(axis.title.x=element_text(size = 12),legend.position="none",
              axis.text.x = element_text(size=12,angle =45, hjust = 1, vjust = 1),
              axis.title.y=element_text(size = 14),axis.text.y = element_text(size=12),
              strip.background  = element_blank(),
              strip.text = element_text(size=12))+
        scale_x_discrete(breaks=unique(coldata$stage),
                         labels=c("","E11.5","","","","E15.5","","","","0dpb","","","","9wpb"))+
        scale_color_manual(values=c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397","#ff6700"))
      ggsave(p1,filename = file.path("~/Pseudo/Result",Species,"Picture",paste0(Num,"Fraction.ExpressedDynamicPseudogene.pdf")),
             device = "pdf",width = 9.5,height = 4.5)
    }
  }
}
