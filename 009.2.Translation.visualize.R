rm(list = ls());gc();rm(list = ls())
if (FALSE) { #has executed
  for (S in c("Human","Mouse")) {
    dir = file.path("~/Pseudo/Data/Seqdata/RPFdb",S)
    Files = grep("metatable$",list.files(dir),value=TRUE)
    filePath <- sapply(Files,function(x){paste(dir,x,sep='/')})   
    data <- lapply(filePath, function(x){ fread(x, header=T,sep = "\t")})  
    a = data[[1]] %>% as.data.frame()
    a %<>% dplyr::filter(.,!grepl("_PAR_Y",Gene_ID)) %>% select(.,-Gene_Name)
    a$Gene_ID %<>% lapply(., function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
    for (i in 2:length(data)) {
      b <- data[[i]] %>% as.data.frame()
      b %<>% dplyr::filter(.,!grepl("_PAR_Y",Gene_ID)) %>% select(.,-Gene_Name)
      b$Gene_ID %<>% lapply(., function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
      a <- merge(a,b,by="Gene_ID")
    }
    rm(data)
    #write.table(a,file.path("~/Pseudo/Data/Seqdata/RPFdb",S,paste0(S,".all.RPKM.txt")),sep = "\t",quote = FALSE)
    #fpkm =a 
    ref <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                    header = FALSE,sep = "\t",stringsAsFactors = F)
    a$type <- ref[match(a$Gene_ID,ref$V4),5]
    a <- dplyr::filter(a,grepl("pseudogene|protein_coding|lncRNA",type,ignore.case = T))
    a$type %<>% gsub(" protein_coding","Coding",.) %>% gsub(" lncRNA","lncRNA",.)
    a[a$type!="Coding" & a$type!="lncRNA","type"] = "Non-dynamic"
    ddg = fread(file.path("~/Pseudo/Result",S,"Savedata/DDG/all.ddg.csv"),header = FALSE)
    a[a$Gene_ID %in% ddg$V1,"type"] = "Dynamic"
    fwrite(a,file.path("~/Pseudo/Data/Seqdata/RPFdb",S,paste0(S,".all.RPKM.txt")),sep = "\t",quote = FALSE,
           row.names = FALSE)
  }
}

for (S in c("Human","Mouse")) {
  a = fread(file.path("~/Pseudo/Data/Seqdata/RPFdb",S,paste0(S,".all.RPKM.txt"))) %>% as.data.frame()
  n = ncol(a)-1
  a$max = apply(a[,2:n],1,max)
  a$median = apply(a[,2:n],1,median)
  a$type %<>% factor(.,levels = c("Coding","Dynamic","Non-dynamic","lncRNA"))
  my_comparisons <- list(c("Coding","Dynamic"),
                         c("Dynamic","Non-dynamic"),
                         c("Non-dynamic","lncRNA"))
  if (S=="Human") {
    m=3
  }else{
    m=4
  }
  p1 = ggplot(a,aes(type,log10(max+1)))+geom_boxplot(aes(fill=type),outlier.colour = "white",notch = T)+
    theme_classic()+ylab("Ribo-seq FPKM (log10)")+
    theme(axis.title.x=element_blank(), axis.title.y=element_text(size=12),
          axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size=10))+
    scale_fill_manual(values = c("#EE1F26","#61439A","#4F69B5","#862461"))+
    stat_compare_means(comparisons = my_comparisons,label.y = c(2.65,2.5,2.35),tip.length = 0.01,size = 3)+
    theme(legend.position='none')+
    coord_cartesian(ylim = c(0,m))
  ggsave(p1,filename = file.path("~/Pseudo/Result/RPFdb",S,"Picture/Max.fpkm.pdf"), 
         device = "pdf",width = 3, height = 4)
  
  tapply(a$median, a$type, mean)
  p2 = ggplot(a,aes(type,median))+geom_boxplot(aes(fill=type),outlier.colour = "white",notch = T)+theme_classic()+ylab("Ribo-seq FPKM (log10)")+theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14),axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1),axis.text.y = element_text(size=12))+
    scale_fill_manual(values = c("#EE1F26","#61439A","#4F69B5","#862461"))+
    stat_compare_means(comparisons = my_comparisons)+
    theme(legend.position='none')+
    coord_cartesian(ylim = c(0,3))
  
  a$num = apply(a[,2:n],1,function(x)sum(x > 1))
  p3 = ggplot(a,aes(type,log10(num+1)))+geom_boxplot(aes(fill=type),outlier.colour = "white",notch = T)+theme_classic()+ylab("Number (FPKM>1, log10)")+theme(axis.title.x=element_blank(), axis.title.y=element_text(size=14),axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1),axis.text.y = element_text(size=12))+
    scale_fill_manual(values = c("#EE1F26","#61439A","#4F69B5","#862461"))+
    stat_compare_means(comparisons = my_comparisons)+
    theme(legend.position='none')
  ggsave(p3,filename = file.path("~/Pseudo/Result/RPFdb",S,"Picture/Num.fpkm1.pdf"), 
         device = "pdf",width = 5, height = 4)
}