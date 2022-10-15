rm(list = ls());gc();rm(list = ls())
S="Human"
wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
a <- fread(file.path(wd,S,"Kaessmann.Mergecount.changed.txt")) %>% as.data.frame()
a[1:3,1:3]
a <- dplyr::filter(a,!grepl("_PAR_Y",V1))
a$V1 <- lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
              header = FALSE,sep = "\t",stringsAsFactors = F)
a$type <- b[match(a$V1,b$V4),5]
a <- dplyr::filter(a,grepl("pseudogene",type,ignore.case = T))
row.names(a) <- a$V1

#Get coldata
if (FALSE) { #has done
  {
    coldata <- data.frame(row.names = colnames(a)[2:(ncol(a)-1)],
                          condition=lapply(colnames(a)[2:(ncol(a)-1)], function(x)strsplit(x,split = "_",fix=TRUE)[[1]][1]) %>% unlist())
    coldata$condition %<>% gsub("KidneyTestis","Kidney",.) %>% 
      gsub("Forebrain","Brain",.) %>% 
      gsub("Hindbrain","Cerebellum",.) %>% factor()
    coldata$stage <- lapply(colnames(a)[2:(ncol(a)-1)], function(x)strsplit(x,split = "_",fix=TRUE)[[1]][2]) %>% unlist()
    
    info <- read.csv("/home/qians/Pseudo/Data/Ref/Human/E-MTAB-6814.sdrf.txt",header = TRUE,sep = "\t",stringsAsFactors = FALSE)
    info2 <- info[,c(1,6:11)]
    info2$Unit..time.unit. <- factor(info2$Unit..time.unit.,levels = c("week","day","month","year"))
    info2 %<>% arrange(Unit..time.unit.,Characteristics.age.)
    info2$time1 <- lapply(info2$Source.Name,function(x)strsplit(x,split = ".",fix=TRUE)[[1]][4]) %>% unlist()
    info2$time2 <- paste0(info2$Characteristics.age.,info2$Unit..time.unit.) %>% gsub("week","wpc",.) %>% gsub("day","dpb",.) %>% gsub("month","mpb",.) %>% gsub("year","ypb",.)
    
    coldata$stage2 <- info2[match(coldata$stage,info2$time1),"time2"]
    coldata$Unit..time.unit. <- info2[match(coldata$stage,info2$time1),"Unit..time.unit."]
    coldata$Unit..time.unit. <- factor(coldata$Unit..time.unit.,levels = c("week","day","month","year"))
    coldata$name <- row.names(coldata)
    coldata %<>% arrange(Unit..time.unit.,stage2)
    coldata$stage3 <- rep(c("S1","S2","S3"),
                          times=c(sum(coldata$Unit..time.unit.=="week"),
                                  sum(coldata$Unit..time.unit.=="day"|coldata$Unit..time.unit.=="month"),
                                  sum(coldata$Unit..time.unit.=="year"))) %>% factor()
    coldata$sex <- lapply(coldata$name, function(x)strsplit(x,split = "_",fix=TRUE)[[1]][3]) %>% unlist() %>% factor() #sample order
    row.names(coldata) <- coldata$name
    write.csv(coldata,file = file.path("~/Pseudo/Result",S,"Savedata","coldata.csv"))
  }
}
coldata = fread(file.path("~/Pseudo/Result",S,"Savedata","coldata.csv"),header = TRUE) %>% as.data.frame()
row.names(coldata) = coldata$V1
coldata <- coldata[,c("condition","stage3","sex")]
#saveRDS(coldata,file = file.path("~/Pseudo/Result",S,"Savedata","dds.rds"))

#expressed pseu
express = fread(file.path("~/Pseudo/Result",S,"Savedata","gene.expressnum.1allfpkm.csv"))  #0.3
express %<>% dplyr::filter(.,num >= 1)
a <- a[express$gene,row.names(coldata)]
all(row.names(coldata)==colnames(a))
dds <- DESeqDataSetFromMatrix(countData = a,
                              colData = coldata,
                              design = ~ condition + stage3 + sex)
dds
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"))
pcaData <- plotPCA(vsd, intgroup=c("condition", "stage3","sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
ggplot(pcaData, aes(PC1, PC2, color=condition,shape=sex)) +
  geom_point(aes(size=stage3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397","#ff6700"))+
  theme_classic()+theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
                        axis.text.x = element_text(size=16),axis.text.y = element_text(size=16),
                        legend.position = "none")
ggsave(filename = file.path("~/Pseudo/Result",S,"Picture","PCA_tissuesexstage.2.pdf"),device = "pdf",
       width = 4,height = 4)
ggplot(pcaData, aes(PC1, PC2, color=condition,shape=sex)) +
  geom_point(aes(size=stage3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c("#3298c8","#33cafe","#c60202","#cb9900","#349a00","#cc3397","#ff6700"))+
  theme_classic()+theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
                        axis.text.x = element_text(size=16),axis.text.y = element_text(size=16),
                        legend.direction = "horizontal",legend.title = element_blank())
ggsave(filename = file.path("~/Pseudo/Result",S,"Picture","PCA_tissuesexstage.pdf"),device = "pdf",
       width = 7.45,height = 4.88)
#save.image(file=file.path("~/Pseudo/Result",S,"Savedata","dds.RData"))