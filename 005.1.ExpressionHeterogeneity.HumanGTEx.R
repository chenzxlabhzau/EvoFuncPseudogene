##bash
#cd ~/Pseudo/Data/Seqdata/GTEx
#wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
#gunzip *gz

#R
############################
##Heatmap
if (TRUE) {
  rm(list = ls())
  Num = "005.1."
  S="Human"
  a <- read.csv("~/Pseudo/Data/Seqdata/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",header = F,sep = "\t",stringsAsFactors = F )
  colnames(a) <- a[3,]
  a <- a[-c(1:3),-c(which(colnames(a)=="Cells - Cultured fibroblasts"),
                    which(colnames(a)=="Cells - EBV-transformed lymphocytes"),
                    which(colnames(a)=="Whole Blood"))]
  a$Name <- lapply(a$Name,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
  b <- lapply(colnames(a),function(x)strsplit(x," ",fixed=T)[[1]][1]) %>% unlist()
  b[which(colnames(a)=="Brain - Cerebellum")] = "Cerebellum"
  b[which(colnames(a)=="Muscle - Skeletal")] = "S. Muscle"
  b[which(colnames(a)=="Minor Salivary Gland")] = "M.Sal.Gland"
  b[which(colnames(a)=="Small Intestine - Terminal Ileum")] = "S. Intestine"
  b[which(colnames(a)=="Nerve - Tibial")] = "NerveTibial"
  b[which(colnames(a)=="Fallopian Tube")] = "FallopianTube"
  b <- b[-c(1:2)]
  c <- t(apply(a[,-c(1:2)],1,function(x){tapply(as.numeric(x),b,mean)})) %>% as.data.frame()
  c$Gene <- a$Name
  c <- apply(c[-ncol(c)],2,function(x){tapply(as.numeric(x),c$Gene,mean)}) %>% as.data.frame() #duplicate gene name
  
  c$Gene = row.names(c)
  d <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  c$type = d[match(c$Gene,d$V4),5]
  c <- dplyr::filter(c,grepl("pseudogene",type))
  
  m1 <- dplyr::filter(c,grepl("polymorphic_pseudogene",type))
  m2 <- dplyr::filter(c,grepl("unprocessed_pseudogene",type))
  m3 <- dplyr::filter(c,grepl("unitary_pseudogene",type))
  m4 <- dplyr::filter(c,Gene %in% setdiff(c$Gene,c(m1$Gene,m2$Gene,m3$Gene)))
  m <- rbind(m1,m2,m3,m4)
  m$type2 = rep(c("Poly","Unp","Uni","Pro"),c(nrow(m1),nrow(m2),nrow(m3),nrow(m4)))
  
  n <- m
  n = apply(n[,-c(ncol(n)-1,ncol(n)-2,ncol(n))], 2,as.numeric) %>% as.data.frame()
  n = apply(n, 2, function(x) 10^6 * x/sum(x)) %>% as.data.frame()
  n$type2 = m$type2
  row.names(n) = m$Gene
  n <- n[apply(n[,-ncol(n)], 1, function(x)sum(x) >0),]
  save(m,n, file = file.path("~/Pseudo/Result/GTEx/Savedata/GTEx.RData"))
  
  row =  n$type2 %>% as.data.frame()
  row.names(row) = row.names(n)
  colnames(row) ="type"
  row$type = factor(row$type,levels = rev(c("Pro","Unp","Poly","Uni")))
  
  n = t(scale(t(n[,-ncol(n)]))) %>% as.data.frame()
  n[n>2]=2
  n[n<=-2]=-2
  ann_colors = list(type=c(Pro="#5da495",Unp="#dfcbe2",Uni="#f7b6ae",Poly="#89c1e7"))
  p=pheatmap(n,annotation_colors = ann_colors,
             annotation_row = row,
             annotation_legend = T,
             border_color = 'black',#设定每个格子边框的颜色，border=F则无边框
             cluster_rows = F, #对行聚类
             cluster_cols = T, #队列聚类
             show_colnames = T, #是否显示列名
             show_rownames = F, #是否显示行名
             main = "Human (GTEx)"
  )
  ggsave(p,filename = file.path("~/Pseudo/Result/GTEx/Picture",paste0(Num,"Heatmap.GTEx.tissue.pdf")),device = "pdf",width = 9,height = 6)
}

##Bar + scatter plot
############################
if (TRUE) {
  rm(list = ls())
  Num = "005.1."
  load(file = "~/Pseudo/Result/GTEx/Savedata/GTEx.RData")
  n <- n[,-ncol(n)]
  n$median <- apply(n[,-ncol(n)],1,median)
  
  c <- c()
  for (i in 1:(ncol(n)-1)) {
    c[i] = sum(apply(n, 1, function(x) x[i] >x[ncol(n)]))
  }
  colnames(n)[which.min(c)]
  b= data.frame(c,colnames(n)[-ncol(n)])
  colnames(b)=c("number","tissue")
  b= b[order(b$number),]
  b$tissue = factor(b$tissue,levels = b$tissue)
  ggplot(b,aes(tissue,number))+geom_bar(stat = "identity",position = "dodge",fill=rep(c("#E39E3E","#C7533B"),c(nrow(b)-1,1)))+theme_classic()+
    ylab("Number of highly expressed pseudogenes")+
    theme(axis.title.x=element_blank(),legend.position="bottom",
          axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1),
          axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
    scale_y_continuous(expand = c(0,0),limits = c(0,6000))
  ggsave(filename = file.path("~/Pseudo/Result/GTEx/Picture",paste0(Num,"Number.HigherexpressedPseu.tissue.pdf")),device = "pdf",width = 8,height = 6)
  
  p1=ggplot(n,aes(log10(median+1),log10(Testis+1)))+geom_point(alpha = 0.1)+
    coord_cartesian(xlim = c(0,4),ylim = c(0,4))+theme_bw()+
    xlab("Mean abundance (All Tissues)")+ylab("Mean abundance (Testis)")+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.title.x=element_text(size=18),
          axis.title.y=element_text(size=18),axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16))+
    geom_abline(slope=1, intercept=0,colour="#6195C8")+
    annotate("text",x=0.55,y=3.8,size=6,colour="#ff6700",
             label=paste(sum(n$Testis > n$median),
                         paste0(round(sum(n$Testis > n$median)*100 / nrow(n)),"%)"),sep = "("))
  p1
  ggsave(p1,filename = file.path("~/Pseudo/Result/GTEx/Picture",paste0(Num,"Scatter.tesis.pdf")),device = "pdf",width = 5,height = 5)
  
  p2=ggplot(n,aes(log10(median+1),log10(Heart+1)))+geom_point(alpha = 0.1)+
    coord_cartesian(xlim = c(0,4),ylim = c(0,4))+theme_bw()+
    xlab("Mean abundance (All Tissues)")+ylab("Mean abundance (Heart)")+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.title.x=element_text(size=18),
          axis.title.y=element_text(size=18),axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16))+
    geom_abline(slope=1, intercept=0,colour="#6195C8")+
    annotate("text",x=3.45,y=0.2,size=6,colour="#c60202",
             label=paste(sum(n$Heart < n$median),
                         paste0(round(sum(n$Heart < n$median)*100 / nrow(n)),"%)"),sep = "("))
  p2
  ggsave(p2,filename = file.path("~/Pseudo/Result/GTEx/Picture","Scatter.heart.pdf"),device = "pdf",width = 5,height = 5)
}

####pseudogene type
############################
if (TRUE) {
  rm(list = ls())
  Num = "005.1."
  library(ggpubr)
  library(patchwork)
  
  load(file = "~/Pseudo/Result/GTEx/Savedata/GTEx.RData")
  n$type2 %<>% factor(.,levels = c("Pro","Unp","Uni","Poly"))
  my_comparisons = list(c("Pro","Uni"),c("Unp","Uni"),c("Poly","Uni"))
  p1 = ggplot(n,aes(type2,log2(Testis+1),fill=type2))+geom_boxplot(notch = TRUE,outlier.alpha = 0)+theme_bw()+
    ylab("Expression level in testis")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))+
    stat_compare_means(comparisons = my_comparisons)+
    scale_fill_manual(values = c("#5da495","#dfcbe2","#f7b6ae","#89c1e7"))+
    guides(fill=FALSE)
  
  n$h80 = apply(n[,-ncol(n)], 1, function(x)sort(x)[24])
  df = data.frame(high= tapply(n$Testis>n$h80, n$type2, sum))
  df$all = table(n$type2) %>% as.numeric()
  df$ratio = df$high/df$all
  df$type = row.names(df)
  p2 = ggplot(df,aes(type,ratio,fill=type))+geom_bar(stat = "identity",position = "dodge")+theme_bw()+
    ylab("Fraction of highly expressed\npseudogenes in testis")+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10))+
    scale_y_continuous(expand = c(0,0),limits = c(0,0.8))+
    scale_fill_manual(values = c("#5da495","#dfcbe2","#f7b6ae","#89c1e7"))+
    guides(fill=FALSE)
    coord_cartesian(ylim = c(0,0))
  ggsave(p1/p2, filename = file.path("~/Pseudo/Result/GTEx/Picture",paste0(Num,"Express.type.testis.pdf")),
         device = "pdf",width = 5,height = 8)
}

####pseudogene age
############################
if (TRUE) {
  rm(list = ls())
  Num = "005.1."
  library(ggpubr)
  library(patchwork)
  species = "Human"
  
  load(file = "~/Pseudo/Result/GTEx/Savedata/GTEx.RData")
  age = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),header = TRUE,sep = ",",stringsAsFactors = FALSE)
 
  n$age = age[match(row.names(n),age[,1]),"age"]
  n %<>% na.omit()
  n$h80 = apply(n[,-ncol(n)], 1, function(x)sort(x)[24])
  df = data.frame(high= tapply(n$Testis>n$h80, n$age, sum))
  df$all = table(n$age) %>% as.numeric()
  df$ratio = df$high/df$all
  df$age = row.names(df) %>% as.numeric()
  
  