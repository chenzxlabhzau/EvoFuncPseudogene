{
  rm(list = ls())
  for (S in c("Human","Mouse")) {
    dir = file.path("/home/qians/Pseudo/Result",S,"Savedata/DDG")
    Files = grep("_ddg.csv$",list.files(dir),value=TRUE)
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
    
    #length
    length <- read.csv(file = file.path("~/Pseudo/Data/Ref",S,paste0(S,".nonreExonlength.bed")),
                       header = FALSE,sep = "\t",stringsAsFactors = FALSE)
    #nonredundant gene length
    gene.length <- tapply(length$V5, length$V4, sum) %>% as.data.frame() 
    gene.length$type = b[match(row.names(gene.length),b$V4),5]
    
    #filter for pseudo-
    gene.length$Gene = row.names(gene.length)
    gene.length %<>% dplyr::filter(.,grepl("pseudogene",type,ignore.case = T))
    
    gene.length$type2 = "Non-dynamic"
    gene.length[gene.length$Gene %in% a$Gene,"type2"] = "Dynamic"
    wilcox.test(gene.length[gene.length$type2=="Dynamic",1],gene.length[gene.length$type2=="Non-dynamic",1],alternative = "greater")
    l =tapply(gene.length$., gene.length$type2, median)
    if (S== "Human") {
      p = ggplot(gene.length,aes(.,fill=type2,color=type2))+geom_density(alpha=0.75,adjust=2)+
        theme_classic()+
        ylab("Density")+xlab("Transcript length (bp)")+
        theme(legend.position = c(0.8, 0.8),legend.title = element_blank(),
              legend.text = element_text(size=8),legend.background = element_blank(),
              legend.key.size = unit(0.4, 'cm'),
              axis.title.x=element_text(size=14),
              axis.text.x = element_text(size=12),
              axis.title.y=element_text(size=14),
              axis.text.y = element_text(size=12))+
        scale_y_continuous(expand = c(0, 0),limits = c(0,1.2))+
        geom_segment(aes(x=l[2], y=1.15, xend=l[1], yend=1.15),color="black")+
        annotate("text", x=mean(l), y=1.16, label="***",size=6)+#annotation_logticks(sides = "b")+
        scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x) ,
          labels = scales::trans_format("log10", scales::math_format(expr = 10^.x)))+
        scale_fill_manual(values = c("#61439A","#4F69B5"))+
        scale_color_manual(values = c("#61439A","#4F69B5"))
    }else{
      p = ggplot(gene.length,aes(.,fill=type2,color=type2))+geom_density(alpha=0.75,adjust=2)+
        theme_classic()+
        ylab("Density")+xlab("Transcript length (bp)")+
        theme(legend.position = c(0.8, 0.8),legend.title = element_blank(),
              legend.text = element_text(size=8),legend.background = element_blank(),
              legend.key.size = unit(0.4, 'cm'),
              axis.title.x=element_text(size=14),
              axis.text.x = element_text(size=12),axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
        scale_y_continuous(expand = c(0, 0))+
        geom_segment(aes(x=l[2], y=1.31, xend=l[1], yend=1.31),color="black")+
        annotate("text", x=mean(l), y=1.33, label="***",size=6)+
        coord_cartesian(xlim = c(10^1,10^4.5),ylim = c(0,1.4))+#annotation_logticks(sides = "b")+
        scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x) ,
      labels = scales::trans_format("log10", scales::math_format(expr = 10^.x)))+
        scale_fill_manual(values = c("#61439A","#4F69B5"))+
        scale_color_manual(values = c("#61439A","#4F69B5"))
      
    }
    ggsave(p,filename = file.path("~/Pseudo/Result",S,"Picture","Density.length.DDG.pdf"),device = "pdf",
          height = 4, width = 3.2)
    
    ##Revise 1
    gene.length$type %<>% gsub("transcribed_","",.) %>% gsub("translated_","",.) %>%
      gsub(" ","",.) %>%
      gsub("unitary_pseudogene","Unitary",.) %>% 
      gsub("polymorphic_pseudogene","Polymorphic",.) %>% 
      gsub("processed_pseudogene","Processed",.)
    gene.length[!gene.length$type %in% c("Unitary","Polymorphic","Processed"),"type"]="Unprocessed"
    table(gene.length$type,gene.length$type2)
    table(gene.length$type,gene.length$type2)[,1]
    table(gene.length$type,gene.length$type2)[,1]/ sum(table(gene.length$type,gene.length$type2)[,1])
    
  }
}
