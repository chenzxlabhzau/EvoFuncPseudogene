rm(list = ls())
library(ChIPseeker)
library(GenomicRanges)
for (S in c("Human","Mouse")) {
  dir = file.path("/home/qians/Pseudo/Data/Ref",S,"TE")
  Files = grep("bed.gz$",list.files(dir),value=TRUE) 
  #Files = Files[-5]
  filePath <- sapply(Files,function(x){paste(dir,x,sep='/')}) %>% as.list()
  
  promoter = fread(file.path("~/Pseudo/Data/Ref",S,paste0(S,".promoter.bed"))) %>% as.data.frame()
  colnames(promoter) = c("chr","start","end","name","strand")
  
  anno = read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                  header = FALSE,sep = "\t",stringsAsFactors = F)
  promoter$type = anno[match(promoter$name,anno$V4),5] %>% gsub(" ","",.)
  promoter %<>% dplyr::filter(.,grepl("pseudogene|protein_coding|lncRNA",type,ignore.case = T))
  promoter$type %<>% gsub("protein_coding","Coding",.)
  promoter[promoter$type!="Coding" & promoter$type!="lncRNA","type"] = "Pseudogene"
  # promoter[promoter$type!="Coding" & promoter$type!="lncRNA","type"] = "Non-dynamic\npseudogene"
  # ddg = fread("~/Pseudo/Result/Human/Savedata/DDG/all.ddg.csv",header = FALSE)
  # promoter[promoter$name %in% ddg$V1,"type"] = "Dynamic\npseudogene"
  
  promoter %<>% na.omit %<>% GRanges()
  promoter = promoter[promoter@ranges@width==3001,]
  
  pro1 = promoter[promoter$type == "lncRNA",]
  pro2 = promoter[promoter$type == "Pseudogene",]
  pro3 = promoter[promoter$type == "Coding",]
  
  tagMatrixList1 <- lapply(filePath, getTagMatrix, windows=pro1)
  tagMatrixList2 <- lapply(filePath, getTagMatrix, windows=pro2)
  tagMatrixList3 <- lapply(filePath, getTagMatrix, windows=pro3)
  
  tagMatrixList = c(tagMatrixList1,tagMatrixList2,tagMatrixList3)
  
  names(tagMatrixList) = rep(c("lncRNA","Pseudogene","Coding"),each=5)
  #plotAvgProf(c(tagMatrixList[n], tagMatrixList[n+5], tagMatrixList[n+10]), xlim=c(-2000, 1000),facet = "col")+theme_classic()
  #plotAvgProf(tagMatrixList, xlim=c(-2000, 1000),facet = "row", conf=0.95,resample=500)
  
  type =lapply(names(tagMatrixList1), function(x)strsplit(x,split = ".",fixed = T)[[1]][3]) %>% unlist()
  type[4] = "TE"
  p = list()
  for (n in 1:5) {
    p[[n]] = plotAvgProf(c(tagMatrixList[n], tagMatrixList[n+5], tagMatrixList[n+10]), xlim=c(-2000, 1000))+theme_classic()+ylab(paste(type[n],"density",sep = " "))+theme(legend.title = element_blank(),axis.title.y=element_text(size=13))
    ggsave(p[[n]],filename = paste0("/home/qians/Pseudo/Result/TEfromUCSC/",S,"/Picture/",
                                    type[n],".pdf"),device = "pdf",width = 5.3,height = 4.5)
  }
}
