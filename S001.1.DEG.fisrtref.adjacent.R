rm(list = ls())
library(data.table)
library(magrittr)
library(stringr)
library(dplyr)
library(DESeq2)

wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
Species = "Human"
a <- fread(file.path(wd,Species,"Kaessmann.Mergecount.changed.txt")) %>% as.data.frame()
coldata <- read.csv(file = file.path("~/Pseudo/Result",Species,"Savedata","coldata.csv"),
                    row.names = 1, sep = ",")
coldata$name %<>% gsub("KidneyTestis","Kidney",.)
coldata$time <- lapply(coldata$stage2,function(x){substr(x,1,str_length(x)-3)}) %>% 
  unlist() %>% as.numeric()

for (i in 1:nrow(coldata)) {
  if (coldata[i,4]=="week") {
    coldata[i,9] = coldata[i,8] * 7
  }
  if (coldata[i,4]=="day") {
    coldata[i,9] = coldata[i,8] * 1 + 280
  }
  if (coldata[i,4]=="month") {
    coldata[i,9] = coldata[i,8] * 30 + 280
  }
  if (coldata[i,4]=="year") {
    coldata[i,9] = coldata[i,8] * 365 + 280
  }
}

coldata$stage4 = paste0(" ",coldata$stage2)
for (i in paste0(" ",c(0,4,6,15,18,19,34,94,127),"dpb")) {
  coldata$stage4 %<>% gsub(i,"Newborn",.,fixed = TRUE)
}
for (i in paste0(" ",c("6mpb","221dpb","226dpb","270dpb","271dpb"))){
  coldata$stage4 %<>% gsub(i,"Infant",.,fixed = TRUE)
}
coldata %>% dplyr::filter(.,grepl("ypb",name)) %>% with(.,table(time))
for (i in paste0(" ",1:4,"ypb")){
  coldata$stage4 %<>% gsub(i,"Toddler",.,fixed = TRUE)
}
for (i in paste0(" ",7:9,"ypb")){
  coldata$stage4 %<>% gsub(i,"SchoolAge",.,fixed = TRUE)
}
for (i in paste0(" ",13:19,"ypb")){
  coldata$stage4 %<>% gsub(i,"Teenager",.,fixed = TRUE)
}
for (i in paste0(" ",25:32,"ypb")){
  coldata$stage4 %<>% gsub(i,"YoungAdult",.,fixed = TRUE)
}
for (i in paste0(" ",39:41,"ypb")){
  coldata$stage4 %<>% gsub(i,"YoungMiddle",.,fixed = TRUE)
}
for (i in paste0(" ",46:55,"ypb")){
  coldata$stage4 %<>% gsub(i,"OldMiddle",.,fixed = TRUE)
}
for (i in paste0(" ",58:63,"ypb")){
  coldata$stage4 %<>% gsub(i,"Senior",.,fixed = TRUE)
}
coldata$stage4 %<>% gsub(" ","",.)

if (all(row.names(coldata) %in% colnames(a)[-1])) {
  a <- dplyr::filter(a,!grepl("_PAR_Y",V1))
  row.names(a) = lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
  a = a[,-1]
  a = a[,row.names(coldata)]
}

deve = read.csv(paste0(wd,Species,"/",Species,".developstage.txt"),
                header = F,sep = "\t",stringsAsFactors = F,row.names = 1)
b = read.csv(paste0("/home/qians/Pseudo/Data/Ref/",Species,"/","gene",Species,".bed"),
             header = F,sep = "\t",stringsAsFactors = FALSE)

for ( Tissue in unique(coldata$condition) ) {
  coldata2 = dplyr::filter(coldata, condition == Tissue & stage4 %in% as.character(deve[Tissue,])) #analysis by tissue
  coldata2$stage4 %<>% factor(.,levels = deve[Tissue,][deve[Tissue,]!=" "] )
  onetissue = a[,row.names(coldata2)] #one tissue matrix
  if (all(row.names(coldata2) == colnames(onetissue))) {
    dds = DESeq2::DESeqDataSetFromMatrix(onetissue,
                                         colData = coldata2,
                                         design = ~stage4)
  }
  dds <- DESeq(dds)
  resultsNames(dds) # lists the coefficients
  #4wpc as ref
  for (s in as.character(deve[Tissue,deve[Tissue,] %>% str_length() >1])[-1]) {
    #c = counts(dds,normalized=T) 
    #all(colnames(c)==row.names(coldata2))
    #mbg = t(apply(c, 1, function(x)tapply(x, coldata2$stage4, mean)))
    res = results(dds,contrast = c("stage4",s,"4wpc")) %>% as.data.frame() %>%
      dplyr::filter(.,log2FoldChange >=0.5 & padj <= 0.05) %>% 
      dplyr::select(.,c(log2FoldChange,padj))
    res$type = b[match(row.names(res),b$V4),5]
    #up
    res %>% dplyr::filter(., grepl("pseudogene",type)) %>% 
      row.names() %>% as.data.frame() %>% 
      write.table(.,file = paste0("/home/qians/Pseudo/Result/",Species,"/Savedata/DEG/4wpcref/",
                                  Species,".",Tissue,".",s,".4wpc.up.pseudogene.txt"),
                  row.names = F,quote = FALSE, col.names = F)
    res %>% dplyr::filter(., grepl("protein_coding",type)) %>% 
      row.names() %>% as.data.frame() %>% 
      write.table(.,file = paste0("/home/qians/Pseudo/Result/",Species,"/Savedata/DEG/4wpcref/",
                                  Species,".",Tissue,".",s,".4wpc.up.coding.txt"),
                  row.names = F,quote = FALSE, col.names = F)
    #down
    res = results(dds,contrast = c("stage4",s,"4wpc")) %>% as.data.frame() %>%
      dplyr::filter(.,log2FoldChange <= -0.5 & padj <= 0.05) %>% 
      dplyr::select(.,c(log2FoldChange,padj))
    res$type = b[match(row.names(res),b$V4),5]
    res %>% dplyr::filter(., grepl("pseudogene",type)) %>% 
      row.names() %>% as.data.frame() %>% 
      write.table(.,file = paste0("/home/qians/Pseudo/Result/",Species,"/Savedata/DEG/4wpcref/",
                                  Species,".",Tissue,".",s,".4wpc.down.pseudogene.txt"),
                  row.names = F,quote = FALSE, col.names = F)
    res %>% dplyr::filter(., grepl("protein_coding",type)) %>% 
      row.names() %>% as.data.frame() %>% 
      write.table(.,file = paste0("/home/qians/Pseudo/Result/",Species,"/Savedata/DEG/4wpcref/",
                                  Species,".",Tissue,".",s,".4wpc.down.coding.txt"),
                  row.names = F,quote = FALSE, col.names = F)
  }
  #adjacent
  for (i in 1:(length(deve[Tissue,deve[Tissue,] %>% str_length() >1])-1)) {
    j= i+1
    si= as.character(deve[Tissue,deve[Tissue,] %>% str_length() >1][i])
    sj= as.character(deve[Tissue,deve[Tissue,] %>% str_length() >1][j])
    
    res = results(dds,contrast = c("stage4",sj,si)) %>% as.data.frame() %>%
      dplyr::filter(.,log2FoldChange >=0.5 & padj <= 0.05) %>% 
      dplyr::select(.,c(log2FoldChange,padj))
    res$type = b[match(row.names(res),b$V4),5]
    #up
    res %>% dplyr::filter(., grepl("pseudogene",type)) %>% 
      row.names() %>% as.data.frame() %>% 
      write.table(.,file = paste0("/home/qians/Pseudo/Result/",Species,"/Savedata/DEG/adjacent/",
                                  Species,".",Tissue,".",sj,".",si,".up.pseudogene.txt"),
                  row.names = F,quote = FALSE, col.names = F)
    res %>% dplyr::filter(., grepl("protein_coding",type)) %>% 
      row.names() %>% as.data.frame() %>% 
      write.table(.,file = paste0("/home/qians/Pseudo/Result/",Species,"/Savedata/DEG/adjacent/",
                                  Species,".",Tissue,".",sj,".",si,".up.coding.txt"),
                  row.names = F,quote = FALSE, col.names = F)
    #down
    res = results(dds,contrast = c("stage4",sj,si)) %>% as.data.frame() %>%
      dplyr::filter(.,log2FoldChange<= -0.5 & padj <= 0.05) %>% 
      dplyr::select(.,c(log2FoldChange,padj))
    res$type = b[match(row.names(res),b$V4),5]
    
    res %>% dplyr::filter(., grepl("pseudogene",type)) %>% 
      row.names() %>% as.data.frame() %>% 
      write.table(.,file = paste0("/home/qians/Pseudo/Result/",Species,"/Savedata/DEG/adjacent/",
                                  Species,".",Tissue,".",sj,".",si,".down.pseudogene.txt"),
                  row.names = F,quote = FALSE, col.names = F)
    res %>% dplyr::filter(., grepl("protein_coding",type)) %>% 
      row.names() %>% as.data.frame() %>% 
      write.table(.,file = paste0("/home/qians/Pseudo/Result/",Species,"/Savedata/DEG/adjacent/",
                                  Species,".",Tissue,".",sj,".",si,".down.coding.txt"),
                  row.names = F,quote = FALSE, col.names = F)
  }
}
