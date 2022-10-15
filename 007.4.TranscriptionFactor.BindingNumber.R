{
  rm(list = ls());gc();rm(list = ls())
  Num = "007.4."
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
    c = read.csv(file.path("/home/qians/Pseudo/Result/GTRDforTF",S,paste0(S,".promoter.number.bed")),header = FALSE,sep = "\t",stringsAsFactors = F)
    c$type= b[match(c$V4,b$V4),5]
    c[grepl("protein_coding",c$type),"type"] = "Coding"
    c[grepl("lncRNA",c$type),"type"] = "lncRNA"
    c[grepl("pseudogene",c$type),"type"] = "Non-dynamic\npseudogene"
    c[c$V4 %in% a$Gene,"type"] = "Dynamic\npseudogene"
    c %<>% dplyr::filter(.,type %in% c("Coding","Non-dynamic\npseudogene","Dynamic\npseudogene","lncRNA"))
    
    random = fread(file.path("~/Pseudo/Result/GTRDforTF/",S,paste0(S,".randompromoter.number.bed"))) %>% as.data.frame()
    random$type = "Intergenic"
    c = rbind(c,random)
    c$type %<>% factor(.,levels = c("Coding","lncRNA","Dynamic\npseudogene","Non-dynamic\npseudogene","Intergenic"))
    #c$type %<>% gsub(" protein_coding","Coding",.) %>% gsub(" lncRNA","lncRNA",.) %>% factor(.,levels = c("Coding","Dynamic\npseudogene","Non-dynamic\npseudogene","lncRNA"))
    my_comparisons <- list(c("Coding","lncRNA"),
                           c("lncRNA","Dynamic\npseudogene"),
                           c("Dynamic\npseudogene","Non-dynamic\npseudogene"),
                           c("Non-dynamic\npseudogene","Intergenic"))
    p = ggplot(c,aes(type,log10(V6+1)))+geom_boxplot(aes(fill=type),outlier.color = "white",notch = T)+
      theme_classic()+ylab("TF enrichment")+
      theme(axis.title.x=element_blank(),legend.position='none',
            axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1),
            axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
      stat_compare_means(comparisons = my_comparisons,label.y = c(3.8,3.6,3.4,3.2),tip.length = 0.015) +
      scale_fill_manual(values = c("#d95f0d","#9ecae1","#61439A","#4F69B5","#862461"))
      #stat_compare_means(label.y = 4) # 添加全局p值
      #stat_compare_means(label = "p.signif",method = "wilcox.test",ref.group = "Non-dynamic\npseudogene",na.rm = TRUE,size=rel(2)) 
    
    ggsave(p,filename = file.path("~/Pseudo/Result/GTRDforTF",S,"Picture",paste0(Num,"Bindsites.pdf")),
           device = "pdf",width = 6,height = 5)
  }
}