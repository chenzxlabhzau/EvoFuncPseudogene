rm(list = ls());gc();rm(list = ls())
Num = "004.4."

# Organ distribution & out of ? organ ----
# * 1.Organ distribution ----
{
  rm(list = ls());gc();rm(list = ls())
  Num = "004.4."
  wd = "/home/qians/Pseudo/Data/Seqdata/Illumina"
  library(dplyr)
  library(tidyr)
  library(ggplotify)
  
  for( species in c("Human","Mouse")) {
    dis.tissue.all = c()
    a <- fread(file.path(wd,species,paste0(species,".allgene.fpkm.txt"))) %>% as.data.frame()
    colnames(a) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.) %>%
      gsub("Forebrain","Brain",.) %>% gsub("Hindbrain","Cerebellum",.)
    b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",species,paste0("gene",species,".bed")), header = FALSE,sep = "\t",stringsAsFactors = F)
    b$V5 %<>% gsub(" ","",.)
    a[,"type"] <- b[match(a$V1,b$V4),5]
    a %<>% dplyr::filter(.,grepl("pseudogene",type,ignore.case = T))
    
    row.names(a) = a$V1
    a = a[,-1]
    a = a[apply(a[,-ncol(a)], 1, max)>0,] #filter 0FPKM in all tissues
    a[1:3,1:3]
    a$sample = colnames(a)[apply(a[,-ncol(a)], 1, which.max)]
    a$tissue = gsub("female_gonad","ovary",a$sample) %>% gsub("male_gonad","testis",.) %>%
      gsub("female_","",.) %>% gsub("male_","",.) %>% 
      lapply(.,function(x)(strsplit(x,"_")[[1]][1])) %>% unlist() %>% Hmisc::capitalize()
    #a$maxexpress = apply(a[,-c(ncol(a),ncol(a)-1,ncol(a)-2)], 1, max)
    dis.tissue = as.data.frame(table(a$tissue)/nrow(a)*100)
    dis.tissue$species = species
    dis.tissue.all = rbind(dis.tissue.all,dis.tissue)
    col = c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397","#ff6700")
    lim = 40
    if (species == "Mouse") {
      col = c("#3298c8","#33cafe","#C2605D","#c60202","#cb9900","#349a00","#cc3397","#ff6700")
      lim = 30
    }
    p1 = ggplot(dis.tissue.all,aes(Var1,Freq,fill=Var1))+geom_bar(stat = "identity",position = "stack")+
      #ggplot(dis.tissue.all,aes(Var2,Freq,fill=Var1))+geom_bar(stat = "identity",position = "fill")+
      theme_classic()+ylab("Tissue distribution (%)")+
      theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1, vjust = 1),
            axis.text.y = element_text(size=10),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=12),
            strip.text.x = element_text(size = 10))+
      scale_y_continuous(expand = c(0, 0),limits = c(0,lim))+
      guides(fill=FALSE)+
      #guides(fill=guide_legend(title=NULL))+
      scale_fill_manual(values = col)
  ggsave(p1, filename = paste0("/home/qians/Pseudo/Result/",species,"/Picture/",Num,"Organdistribution.pdf"),
         device = "pdf",width = 2.5, height = 3)
  }
}

# * Organ distribution (age) ----
if (FALSE) {
  {
    rm(list = ls())
    Num = "004.4."
    wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
    library(dplyr)
    library(tidyr)
    library(ggplotify)
    dis.tissue.all = c()
    for( species in c("Human","Mouse")) {
      a <- fread(paste0(wd,species,"/",species,".allgene.fpkm.txt")) %>% as.data.frame()
      colnames(a) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.) %>%
        gsub("Forebrain","Brain",.) %>% gsub("Hindbrain","Cerebellum",.)
      b <- read.csv(file.path("~/MamDC/Data/Ref",species,paste0("gene",species,".bed")),header = F,sep = "\t") 
      b$V5 %<>% gsub(" ","",.)
      a[,"type"] <- b[match(a$V1,b$V4),5]
      a %<>% dplyr::filter(.,grepl("pseudogene",type,ignore.case = T))
      
      row.names(a) = a$V1
      a = a[,-1]
      a = a[apply(a[,-ncol(a)], 1, max)>0,] #filter 0FPKM in all tissues
      a[1:3,1:3]
      a$sample = colnames(a)[apply(a[,-ncol(a)], 1, which.max)]
      a$tissue = lapply(a$sample,function(x)(strsplit(x,"_")[[1]][1])) %>% unlist()
      
      age = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),header = TRUE,sep = ",",stringsAsFactors = FALSE)
      a$age = age[match(row.names(a),age[,1]),"age"]
      a %<>% na.omit() %>% dplyr::filter(.,age > -11)
      #a$group = cut(a$age,breaks = c(-450,-300,-90,0),labels =  c("Old","Middle","Young"))
      
      
      dis.tissue = as.data.frame(table(a$tissue)/nrow(a))
      dis.tissue$species = species
      dis.tissue.all = rbind(dis.tissue.all,dis.tissue)
    }
    p1 = ggplot(dis.tissue.all,aes(Var1,Freq,fill=Var1))+geom_bar(stat = "identity",position = "stack")+
      #ggplot(dis.tissue.all,aes(Var2,Freq,fill=Var1))+geom_bar(stat = "identity",position = "fill")+
      theme_classic()+facet_grid(.~species)+ylab("Tissue distribution")+
      theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1, vjust = 1),
            axis.text.y = element_text(size=10),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=12),
            strip.text.x = element_text(size = 10))+
      guides(fill=FALSE)+
      #guides(fill=guide_legend(title=NULL))+
      scale_fill_manual(values = c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397","#ff6700"))
    ggsave(p1, filename = paste0("/home/qians/Pseudo/Result/Human/Picture/",Num,"HM.Organdistribution.pdf"),
           device = "pdf",width = 4, height = 3)
  }
}

# * Out of organ ----
{
  rm(list = ls())
  Num = "004.4."
  wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
  library(dplyr)
  library(tidyr)
  library(ggplotify)
  for( species in c("Human","Mouse")) {
    a <- fread(paste0(wd,species,"/",species,".allgene.fpkm.txt")) %>% as.data.frame()
    colnames(a) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.) %>%
      gsub("Forebrain","Brain",.) %>% gsub("Hindbrain","Cerebellum",.)
    b <- read.csv(file.path("~/MamDC/Data/Ref",species,paste0("gene",species,".bed")),header = F,sep = "\t") 
    b$V5 %<>% gsub(" ","",.)
    a[,"type"] <- b[match(a$V1,b$V4),5]
    a %<>% dplyr::filter(.,grepl("pseudogene",type,ignore.case = T))
    
    row.names(a) = a$V1
    a = a[,-1]
    a = a[apply(a[,-ncol(a)], 1, max)>0,] #filter 0FPKM in all tissues
    a[1:3,1:3]
    group = lapply(colnames(a)[-ncol(a)], function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][1]) %>% unlist() %>%
      gsub("Cerebellum","Brain",.)
    mbg = apply(a[,-ncol(a)], 1, function(x)tapply(x, group, max)) %>% t %>% as.data.frame()
    mbg$sum = apply(mbg,1,sum)
    mbg = apply(mbg[,1:(ncol(mbg)-1)], 2, function(x)x/mbg$sum) %>% as.data.frame() #tissue specificity
    
    age = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),header = TRUE,sep = ",",stringsAsFactors = FALSE)
    mbg$age = age[match(row.names(mbg),age[,1]),"age"]
    mbg %<>% na.omit()
    
    pair = read.csv("~/Pseudo/Data/Ref/Human/pseudo.coding.pair.txt",header = T,sep = "\t")
    mbg = mbg[pair$Pseudogene,] %>% na.omit()
    
    ta = apply(mbg[,-ncol(mbg)], 2,function(x)tapply(x, mbg$age, mean)) %>% as.data.frame()
    ta$age = row.names(ta) %>% as.numeric()
    ta %<>% pivot_longer(.,cols=1:(ncol(mbg)-1))
    ggplot(ta,aes(age,value,group=name,color=name))+geom_point()+geom_smooth(method = "lm")+
      scale_color_manual(values = c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397","#ff6700"))
    
    apply(mbg, 2,function(x)cor.test(x,mbg$age,method = "spearman"))
    
    a$sample = colnames(a)[apply(a[,-ncol(a)], 1, which.max)]
    a$tissue = lapply(a$sample,function(x)(strsplit(x,"_")[[1]][1])) %>% unlist()
    #a$maxexpress = apply(a[,-c(ncol(a),ncol(a)-1,ncol(a)-2)], 1, max)
    dis.tissue = as.data.frame(table(a$tissue)/nrow(a))
    dis.tissue$species = species
    dis.tissue.all = rbind(dis.tissue.all,dis.tissue)
  }
  p1 = ggplot(dis.tissue.all,aes(Var1,Freq,fill=Var1))+geom_bar(stat = "identity",position = "stack")+
    #ggplot(dis.tissue.all,aes(Var2,Freq,fill=Var1))+geom_bar(stat = "identity",position = "fill")+
    theme_classic()+facet_grid(.~species)+ylab("Tissue distribution")+
    theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1, vjust = 1),
          axis.text.y = element_text(size=10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=12),
          strip.text.x = element_text(size = 10))+
    guides(fill=FALSE)+
    #guides(fill=guide_legend(title=NULL))+
    scale_fill_manual(values = c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397","#ff6700"))
  ggsave(p1, filename = paste0("/home/qians/Pseudo/Result/Human/Picture/",Num,"HM.Organdistribution.pdf"),
         device = "pdf",width = 4, height = 3)
}


age = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),header = TRUE,sep = ",",stringsAsFactors = FALSE)
a$age = age[match(row.names(a),age[,1]),"age"]
ta = na.omit(a)