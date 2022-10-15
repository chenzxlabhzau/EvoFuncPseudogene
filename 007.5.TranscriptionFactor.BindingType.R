{
  rm(list = ls());gc();rm(list = ls())
  Num = "007.5."
  for (S in c("Human","Mouse")) {
    dir = file.path("/home/qians/Pseudo/Result",S,"Savedata/DDG")
    Files = grep("ddg.csv$",list.files(dir),value=TRUE)
    filePath <- sapply(Files,function(x){paste(dir,x,sep='/')})   
    data <- lapply(filePath, function(x){ read.csv(x, header=TRUE,sep = ",",stringsAsFactors = FALSE)})  
    a <- data[1] %>% as.data.frame()
    a = a[,1,drop=FALSE]
    colnames(a) <- c("Gene")
    for (i in 2:length(data)) {
      b <- data[i] %>% as.data.frame()
      colnames(b) <- c("Gene")
      b= b[,1,drop=FALSE]
      a <- rbind(a,b)
    }
    b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                  header = FALSE,sep = "\t",stringsAsFactors = F)
    c =fread(file.path("~/Pseudo/Result/GTRDforTF",S,paste0(S,".promoter.type3.bed")),header = FALSE)
    
    c$type= b[match(c$V2,b$V4),5]
    c[grepl("protein_coding",c$type),"type"] = "Coding"
    c[grepl("lncRNA",c$type),"type"] = "lncRNA"
    c[grepl("pseudogene",c$type),"type"] = "Non-dynamic\npseudogene"
    c[c$V2 %in% a$Gene,"type"] = "Dynamic\npseudogene"
    c %<>% dplyr::filter(.,type %in% c("Coding","Non-dynamic\npseudogene","Dynamic\npseudogene","lncRNA"))
    c$type %<>% factor(.,levels = c("Coding","lncRNA","Dynamic\npseudogene","Non-dynamic\npseudogene"))
    my_comparisons <- list(c("Coding","lncRNA"),
                           c("lncRNA","Dynamic\npseudogene"),
                           c("Dynamic\npseudogene","Non-dynamic\npseudogene"))
    p=ggplot(c,aes(type,log10(V1+1)))+geom_boxplot(aes(fill=type),outlier.color = "white",notch = T)+
      theme_classic()+ylab("TF diversity")+
      theme(axis.title.x=element_blank(),legend.position='none',
            axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1),
            axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
      stat_compare_means(comparisons = my_comparisons,label.y = c(3.1,3,2.9),tip.length = 0.015) +
      scale_fill_manual(values = c("#d95f0d","#9ecae1","#61439A","#4F69B5"))
    ggsave(p,filename = file.path("~/Pseudo/Result/GTRDforTF",S,"Picture",paste0(Num,"Bindtype.pdf")),
           device = "pdf",width = 6,height = 5)
  }
}
