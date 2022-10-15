
# Tissue specificity ------------------------------------------------------
{
  rm(list = ls());gc();rm(list = ls())
  wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
  for( Species in c("Human","Mouse")) {
    #for (Sex in c("Female","Male")) {
    a <- fread(file.path(wd,Species,paste0(Species,".allgene.fpkm.txt"))) %>% as.data.frame()
    if (Species == "Mouse") {
      a <- fread(file.path(wd,Species,paste0(Species,".allgene.fpkm.Kaessmann.txt"))) %>% as.data.frame()
    }
    a[1:3,1:3]
    #a <- dplyr::filter(a,!grepl("_PAR_Y",V1))
    #a$V1 <- lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
    b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",Species,paste0("gene",Species,".bed")), header = FALSE,sep = "\t",stringsAsFactors = F)
    a$type <- b[match(a$V1,b$V4),5]
    a <- dplyr::filter(a,grepl("pseudogene|protein_coding|lncRNA",type,ignore.case = T))
    row.names(a) <- a$V1
    a <- a[,-c(1,ncol(a))]
    a[1:3,1:3]
    #a <- apply(a, 2, function(x){10^6 * x/sum(x)}) %>% as.data.frame()
    colnames(a) <- gsub("KidneyTestis_4w_Male","Kidney_4w_Male",colnames(a))
    coldata <- read.csv(file = file.path("~/Pseudo/Result",Species,"Savedata","coldata.csv"),
                        row.names = 1, sep = ",")
    
    g <- c()
    for (Tissue in unique(sort(coldata$condition))) {
      #for (Sex in c("Female","Male")) {
      c <- coldata[coldata$condition==Tissue,]
      row.names(c) <- gsub("KidneyTestis","Kidney",row.names(c))
      d <- a[,row.names(c)]
      #d <- dplyr::select(a,starts_with(Tissue)) # %>% dplyr::select(.,ends_with(Sex,ignore.case = FALSE))
      c$ts <- paste(c$condition,c$stage2,sep = "_")
      d <- t(apply(d,1,function(x){tapply(x,c$ts,max)})) %>% as.data.frame()
      tau <- function(x){
        x <- x[apply(x, 1, function(y)sum(y > 0.3) >= 1),]
        x <-  as.data.frame(t(apply(x, 1,function(y){1-y/max(y)})))
        x$tau <- apply(x,1,function(y)sum(y)/(length(y)-1))
        return(as.data.frame(x))
      }
      d <- tau(d)
      d$tissue <- Tissue
      d$type <- b[match(row.names(d),b$V4),5]
      g <- rbind(g,d[,c("type","tau","tissue")])
    }
    g$type <- lapply(g$type,function(x)strsplit(x," ")[[1]][2])  %>% unlist()
    g[g$type=="protein_coding","type"]="Coding"
    g[grepl("pseudogene",g$type),"type"]="Pseudogene"
    g$type <- factor(g$type,levels = c("Coding","Pseudogene","lncRNA"))
    
    testP1 = c()
    for (i in unique(g$tissue)) {
      test = wilcox.test(g[g$tissue==i&g$type=="Coding","tau"],g[g$tissue==i&g$type=="Pseudogene","tau"])
      testP1 = c(testP1,test$p.value)
    }
    
    testP2 = c()
    for (i in unique(g$tissue)) {
      test = wilcox.test(g[g$tissue==i&g$type=="lncRNA","tau"],g[g$tissue==i&g$type=="Pseudogene","tau"])
      testP2 = c(testP2,test$p.value)
    }
    
    ggplot(g,aes(type,tau))+geom_boxplot(aes(fill=type),notch = TRUE,outlier.colour = "white")+
      ylab("Stage specificity")+theme_classic()+ 
      theme(axis.title.x=element_blank(),legend.position='none',
            axis.text.x = element_text(size=16,angle = 45, hjust = 1, vjust = 1),
            axis.title.y=element_text(size=18),axis.text.y = element_text(size=16))+
      scale_fill_manual(values = c("#d95f0d","#fc9272","#9ecae1"))+ 
      facet_grid(. ~ tissue )+
      geom_text(aes(x,y,label =lab),
                data = data.frame(x = 1,
                                  y = -0.03,
                                  lab = cut(testP1,breaks = c(1,0.05,0.01,0.001,-0.1),
                                            labels = rev(c("n.s.","*","**","***"))),
                                  tissue = factor(unique(g$tissue)),levels = unique(g$tissue)),
                vjust = 1, size=5)+
      geom_text(aes(x,y,label =lab),
                data = data.frame(x = 3,
                                  y = ((testP2>=0.05)/50)-0.03,
                                  lab = cut(testP2,breaks = c(1,0.05,0.01,0.001,-0.1),
                                            labels = rev(c("n.s.","*","**","***"))),
                                  tissue = factor(unique(g$tissue)),levels = unique(g$tissue)),
                vjust = 1, size=5)
    
    ggsave(filename = file.path("~/Pseudo/Result",Species,"Picture",paste0(Species,".stage.specificity.pdf")),
           device = "pdf",width = 8.2,height = 5)
    
    tmp = g %>% dplyr::filter(.,type %in% "Pseudogene")
    tmp$type2 = b[match(row.names(tmp),b$V4),5]
    tmp[grepl("unprocessed_pseudogene",tmp$type2),"type2"] = "Unprocessed"
    tmp[grepl("polymorphic_pseudogene",tmp$type2),"type2"] = "Polymorphic"
    tmp[grepl("unitary_pseudogene",tmp$type2),"type2"] = "Unitary"
    tmp[!grepl("Unprocessed|Polymorphic|Unitary",tmp$type2),"type2"] = "Processed"
    ggplot(tmp,aes(type2,tau))+geom_boxplot()+facet_grid(.~tissue)
  }
}


