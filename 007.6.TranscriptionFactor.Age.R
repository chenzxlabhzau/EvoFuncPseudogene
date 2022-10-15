rm(list = ls());gc();rm(list = ls())
Num = "007.6"

# TF number -------------------
for (species in c("Human","Mouse")) {
  a = read.csv(file.path("/home/qians/Pseudo/Result/GTRDforTF",species,paste0(species,".promoter.number.bed")),
               header = FALSE,sep = "\t",stringsAsFactors = F)
  age = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),
                 header = TRUE,sep = ",",stringsAsFactors = FALSE)
  if (species=="Human") {
    age$number = a[match(age[,1],a$V4),5]
  }else {
    age$number = a[match(age[,1],a$V4),6]
  }
  cor = cor.test(age$age,age$number,method = "spearman")
  age$age %<>% as.factor() %>% as.numeric()



  
  p = ggplot(age,aes(age,log10(number+1),group=age,color=age))+geom_boxplot(outlier.colour = "white",notch = TRUE)+
    geom_smooth(se=TRUE, aes(group=1),method = "lm")+
    scale_color_gradient(high = "#D6B6BB",low = "#9c0b1f")+
    theme_bw()+xlab("From old to young")+ylab("log10(TF-binding sites)")+
    theme(axis.text.x = element_blank(),axis.title.x = element_text(size=16),
          axis.text.y = element_text(size=14),axis.title.y = element_text(size=16),
          axis.ticks.x = element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
    guides(color=FALSE)+
    annotate(geom="text", x=3, y=0.1, size=6,
             label=paste0("rho=",round(cor$estimate,2),",P=",round(cor$p.value,4)))
  ggsave(p, filename = paste0("~/Pseudo/Result/GTRDforTF/",species,"/Picture/",Num,".",species,".Bindsites.age.pdf"),
         width = 5,height = 4)
  
  #add type
  b = read.csv(file.path("/home/qians/Pseudo/Data/Ref",species,paste0("gene",species,".bed")),
               header = FALSE,sep = "\t",stringsAsFactors = F)
  m1 <- dplyr::filter(b,grepl("polymorphic_pseudogene",V5))
  m2 <- dplyr::filter(b,grepl("unprocessed_pseudogene",V5))
  m3 <- dplyr::filter(b,grepl("unitary_pseudogene",V5))
  m4 <- dplyr::filter(b,V4 %in% setdiff(b$V4,c(m1$Gene,m2$Gene,m3$Gene)))
  m <- rbind(m1,m2,m3,m4)
  m$type2 = rep(c("Polymorphic","Unprocessed","Unitary","Processed"),c(nrow(m1),nrow(m2),nrow(m3),nrow(m4)))
  age$type = m[match(age[,1],m$V4),"type2"]
  
  options(digits = 2)
  df = matrix(0,nrow = 4,ncol = 2) %>% as.data.frame()
  for (i in 1:4) {
    cor = with(age[age$type==unique(age$type)[i],],cor.test(number,as.numeric(age),method = "spearman"))
    df[i,1] = cor$estimate
    df[i,2] = cor$p.value
  }
  row.names(df) = unique(age$type)
  df$V3 = paste0("rho=",round(df$V1,2),", P=",format(df$V2,2))
  age$age %<>% as.factor()
  p = ggplot(age,aes(age,log10(number+1)))+geom_boxplot(outlier.colour = "white",notch = TRUE)+
    geom_smooth(se=TRUE, aes(group=1),method = "lm")+
    scale_color_gradient(high = "#D6B6BB",low = "#9c0b1f")+
    theme_bw()+xlab("From old to young")+ylab("Number of state")+
    theme(axis.text.x = element_blank(),axis.title.x = element_text(size=16),
          axis.text.y = element_text(size=14),axis.title.y = element_text(size=16),
          axis.ticks.x = element_blank(),
          panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          strip.text.x = element_text(size = 12),strip.background = element_rect(fill="white"))+
    guides(color=FALSE)+
    facet_wrap(.~type,nrow = 2)+
    geom_text(aes(x,y,label =lab),
              data = data.frame(x = 3,
                                y = 0.3,
                                lab = as.character(df$V3),
                                type = factor(unique(age$type)),levels = unique(age$type)), #facet use name
              vjust = 1, size=4)
  ggsave(p, filename = paste0("~/Pseudo/Result/GTRDforTF/",species,"/Picture/",Num,".",species,".Bindsites.age.type.pdf"),
         width = 7,height = 6)
}  

# TF type -------------------
##* binding type (incomplete) ----------------
rm(list = ls())
for (species in c("Human","Mouse")) {
  a = fread(file.path("~/Pseudo/Result/GTRDforTF",species,paste0(species,".promoter.type3.bed")),
            header = FALSE) %>% as.data.frame()
  #Human.promoter.type3.bed
  age = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),
                 header = TRUE,sep = ",",stringsAsFactors = FALSE)
  if (species=="Human") {
    age$number = a[match(age[,1],a$V2),1]
  }else {
    age$number = a[match(age[,1],a$V4),6]
  }
  cor = cor.test(age$age,age$number,method = "spearman")
  age$age %<>% as.factor() %>% as.numeric()
}