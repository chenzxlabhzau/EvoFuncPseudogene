
# Kaessmann ---------------------------------------------------------------
{
  rm(list = ls());gc();rm(list = ls())
  S="Mouse"
  wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
  a <- fread(file.path(wd,S,"Kaessmann.Mergecount.changed.txt")) %>% as.data.frame()
  a[1:3,1:3]
  # a <- dplyr::filter(a,!grepl("_PAR_Y",V1))
  #a$V1 <- lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist() #delete version
  row.names(a) = a$V1
  a = a[,-1]
  convert = fread("/NAS/qians/geo_submission/mouse_reference_iso.txt",
                  header = FALSE) %>% as.data.frame()
  a1 = a[grepl("ENSMUSG",row.names(a)),]
  a2 = a[grepl("PB",row.names(a)),]
  a2$name = convert[match(row.names(a2),convert$V2),1]
  a2 %<>% na.omit()
  mbg = apply(a2[,-ncol(a2)], 2, function(x)tapply(x, a2$name, max)) %>% as.data.frame()
  a = rbind(a1,mbg)
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  a$type <- b[match(row.names(a),b$V4),5]
  a <- dplyr::filter(a,grepl("pseudogene",type,ignore.case = T))
  
  if (FALSE) { # done
    coldata <- data.frame(row.names = colnames(a)[1:(ncol(a)-1)],
                          condition=lapply(colnames(a)[1:(ncol(a)-1)], function(x)strsplit(x,split = "_",fix=TRUE)[[1]][1]) %>% unlist())
    coldata$stage <- lapply(colnames(a)[2:(ncol(a)-1)],function(x)strsplit(x,split = "_")[[1]][2]) %>%
      unlist() %>% sub("s","",.,fixed = TRUE) #stage has .5
    for (i in 1:nrow(coldata)) {
      if(coldata[i,"stage"] %in% as.character(10:18)){coldata[i,"stage"] = paste0("e",coldata[i,"stage"],".5")}
    }
    coldata$stage <- factor(coldata$stage,levels = c("e10.5","e11.5","e12.5","e13.5","e14.5","e15.5","e16.5",
                                                     "e17.5","e18.5","0dpb","3dpb","2wpb","4wpb","9wpb"))
    coldata$name <- row.names(coldata)
    coldata %<>% arrange(coldata$stage) #note order
    coldata$sex <- lapply(coldata$name,function(x)strsplit(x,split = "_")[[1]][3]) %>%
      unlist() %>% factor()
    row.names(coldata) <- coldata$name
    coldata$stage2 <- rep(c("S1","S2","S3"),
                          times=c(sum(coldata$stage %in% c("e10.5","e11.5","e12.5","e13.5","e14.5",
                                                           "e15.5","e16.5","e17.5","e18.5")),
                                  sum(coldata$stage %in% c("0dpb","3dpb","2wpb")),
                                  sum(coldata$stage %in% c("4wpb","9wpb")))) %>% factor()
    coldata$condition %<>% factor()        
    write.csv(coldata,file = file.path("~/Pseudo/Result",S,"Savedata","coldata.csv"),quote = F)
  }
  
  coldata = fread( file.path("~/Pseudo/Result",S,"Savedata","coldata.csv"),header = TRUE) %>% as.data.frame()
  row.names(coldata) = coldata$V1
  coldata <- coldata[,c("condition","stage2","sex")]
  
  #expressed pseu
  express = fread(file.path("~/Pseudo/Result",S,"Savedata","gene.expressnum.1allfpkm.csv")) 
  express %<>% dplyr::filter(.,num >= 1)
  a <- a[express$gene,row.names(coldata)]
  
  #count matrix and colnum data should be in same order
  if (all(rownames(coldata) == colnames(a))) {
    dds <- DESeqDataSetFromMatrix(countData = a,
                                  colData = coldata,
                                  design = ~ condition + stage2 + sex)
    dds
    dds <- DESeq(dds)
    #vsd <- vst(dds, blind=FALSE)
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    #saveRDS(c(coldata,vsd),file = file.path("~/Pseudo/Result",S,"Savedata","dds.rds"))
    plotPCA(vsd, intgroup=c("condition"))
    pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
    pcaData <- plotPCA(vsd, intgroup=c("condition", "stage2","sex"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes(PC1, PC2, color=condition,shape=sex)) +
      geom_point(aes(size=stage2)) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      scale_color_manual(values = c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397","#ff6700"))+
      theme_classic()+theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
                            axis.text.x = element_text(size=16),axis.text.y = element_text(size=16),
                            legend.direction = "horizontal",legend.title = element_blank())
    ggsave(filename = file.path("~/Pseudo/Result",S,"Picture","PCA_tissuesexstage.Kaessmann.pdf"),
           device = "pdf",width = 7,height = 4)
  }
}


# Our data ---------------------------------------------------------------
{
  rm(list = ls());gc();rm(list = ls())
  S="Mouse"
  wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
  a <- fread(file.path(wd,S,"Ourdata.Mergecount.txt")) %>% as.data.frame()
  a[1:3,1:3]
  row.names(a) = a$V1
  a = a[,-1]
  
  convert = fread("/NAS/qians/geo_submission/mouse_reference_iso.txt",
                  header = FALSE) %>% as.data.frame()
  a1 = a[grepl("ENSMUSG",row.names(a)),]
  a2 = a[grepl("PB",row.names(a)),]
  a2$name = convert[match(row.names(a2),convert$V2),1]
  a2 %<>% na.omit()
  mbg = apply(a2[,-ncol(a2)], 2, function(x)tapply(x, a2$name, max)) %>% as.data.frame()
  a = rbind(a1,mbg)
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  a$type <- b[match(row.names(a),b$V4),5]
  a <- dplyr::filter(a,grepl("pseudogene",type,ignore.case = T))
  
  coldata2 = colnames(a)[-ncol(a)] %>% as.data.frame()
  coldata2$condition = coldata2$. %>% gsub("female_gonad","female_ovary",.) %>% 
    gsub("male_gonad","male_testis",.) %>% 
    lapply(., function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][2]) %>% unlist() %>% capitalize()
  coldata2$stage2 = "S3"
  coldata2$sex = coldata2$. %>% lapply(., function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][1]) %>% 
    unlist() %>% capitalize()
  row.names(coldata2) = coldata2$.
  coldata2 = coldata2[,-1]
  
  #expressed pseu
  express = fread(file.path("~/Pseudo/Result",S,"Savedata","gene.expressnum.1allfpkm.csv")) 
  express %<>% dplyr::filter(.,num >= 1)
  a <- a[express$gene,row.names(coldata2)]
  
  #count matrix and colnum data should be in same order
  if (all(rownames(coldata2) == colnames(a))) {
    dds <- DESeqDataSetFromMatrix(countData = a,
                                  colData = coldata2,
                                  design = ~ condition + sex)
    dds
    dds <- DESeq(dds)
    vsd <- vst(dds, blind=FALSE)
    #vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    #saveRDS(c(coldata,vsd),file = file.path("~/Pseudo/Result",S,"Savedata","dds.rds"))
    plotPCA(vsd, intgroup=c("condition"))
    pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
    pcaData <- plotPCA(vsd, intgroup=c("condition", "sex"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes(PC1, PC2, color=condition,shape=sex)) +
      geom_point(size=5) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      scale_color_manual(values = c("#3298c8","#33cafe","#349a00","#c60202","#cc3397","#ff6700"))+
      theme_classic()+theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
                            axis.text.x = element_text(size=16),axis.text.y = element_text(size=16),
                            legend.direction = "horizontal",legend.title = element_blank())
    ggsave(filename = file.path("~/Pseudo/Result",S,"Picture","PCA_tissuesexstage.Our.pdf"),
           device = "pdf",width = 7,height = 4)
  }
}