rm(list = ls());gc();rm(list = ls())
Num = "004.1."

library(ggpubr)
for( Species in c("Human","Mouse")) {
  #for (Sex in c("Female","Male")) {
  wd = "/home/qians/Pseudo/Data/Seqdata/Illumina"
  a <- fread(file.path(wd,Species,paste0(Species,".allgene.fpkm.txt"))) %>% as.data.frame()
  a[1:3,1:3]
  #a <- dplyr::filter(a,!grepl("_PAR_Y",V1))
  #a$V1 <- lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",Species,paste0("gene",Species,".bed")), header = FALSE,sep = "\t",stringsAsFactors = F)
  a$type <- b[match(a$V1,b$V4),5]
  a <- dplyr::filter(a,grepl("pseudogene|protein_coding|lncRNA",type,ignore.case = T))
  row.names(a) <- a$V1
  a <- a[,-c(1,ncol(a))]
  a[1:3,1:3]
  #a <- apply(a, 2, function(x){10^6 * x/sum(x)})
  colnames(a) <- gsub("KidneyTestis_4w_Male","Kidney_4w_Male",colnames(a))
  sample <- colnames(a) %>% as.data.frame()
  sample$tissue <- lapply(as.character(sample$.),function(x)(strsplit(x,"_")[[1]][1])) %>%  
    gsub("KidneyTestis","Kidney",.) %>% gsub("Forebrain","Brain",.) %>% 
    gsub("Hindbrain","Cerebellum",.) %>% unlist()
  sample$sex <- lapply(as.character(sample$.),function(x)(strsplit(x,"_")[[1]][3])) %>% unlist()
  sample$TS <- paste(sample$tissue,sample$sex,sep = "_")
  
  all(colnames(a)==sample$.)
  b1 <- t(apply(a,1,function(x){tapply(x,sample$tissue,max)})) %>% as.data.frame() 
  #tissue or tissue + sex
  #b1 %<>% dplyr::select(.,ends_with(Sex,ignore.case = FALSE))
  tau <- function(x){
    x <- x[apply(x, 1, function(y)sum(y >1 ) >=1),]
    x <-  as.data.frame(t(apply(x, 1,function(y){1-y/max(y)})))
    x$tau <- apply(x,1,function(y)sum(y)/(length(y)-1))
    return(as.data.frame(x))
  }
  b1 <- tau(b1)
  b1$type <- b[match(row.names(b1),b$V4),5]
  b1$type <- lapply(b1$type,function(x)strsplit(x," ")[[1]][2])  %>% unlist()
  ######################
  #Tissue specificity ~ Coding, Pseudogene and lncRNA
  b1[b1$type=="protein_coding","type"]="Coding"
  b1[grepl("pseudogene",b1$type),"type"]="Pseudogene"
  b1$type <- factor(b1$type,levels = c("Coding","Pseudogene","lncRNA"))
  #my_comparisons = list(c("Coding","Pseudogene"),c("Pseudogene","lncRNA"))
  ggplot(b1,aes(type,tau,fill=type))+geom_boxplot(notch = T,outlier.colour = "white")+theme_classic()+ylab("Tissue specificity")+
    theme(axis.title.x=element_blank(), 
          axis.title.y=element_text(size=18),
          axis.text.x = element_text(size=16,angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size=16))+
    #scale_fill_manual(values = c("#EE1F26","#862461","#61439A","#4F69B5"))+ 
    #scale_fill_manual(values = c("#C1C2E2","#8CB78D","#525252"))+ 
    scale_fill_manual(values =c("#d95f0d","#fc9272","#9ecae1")) +
    #stat_compare_means(comparisons = my_comparisons)+
    geom_segment(aes(x=1.05, y=1, xend=1.95, yend=1),size=0.2)+
    annotate("text", x=1.5, y=1.01, label="***",size=6)+
    geom_segment(aes(x=2.05, y=1, xend=2.95, yend=1),size=0.2)+
    annotate("text", x=2.5, y=1.01, label="***",size=6)+
    theme(legend.position='none')
  ggsave(filename = file.path("~/Pseudo/Result",Species,"Picture",paste0(Num,Species,".tissue.specificity.pdf")),device = "pdf",
         width = 4.5,height = 4.5)
  ######################
}

