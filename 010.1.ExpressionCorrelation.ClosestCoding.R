#2020.08.04
#########################################
#bash
cd ~/Pseudo/Data/Ref/Human
grep protein_coding geneHuman.bed > geneHuman.coding.bed
grep pseudogene geneHuman.bed > geneHuman.pseudogene.bed
#bedtools v2.25.0
bedtools closest -a geneHuman.pseudogene.bed -b geneHuman.coding.bed -io -t first -D b > pseudo.coding.closest.bed

cd ~/Pseudo/Data/Ref/Mouse
grep protein_coding geneMouse.bed > geneMouse.coding.bed
grep pseudogene geneMouse.bed > geneMouse.pseudogene.bed
bedtools closest -a geneMouse.pseudogene.bed -b geneMouse.coding.bed -io -t first -D b > pseudo.coding.closest.bed
#########################################

#########################################
#R
##Distance
rm(list = ls())
wd = "/home/qians/Pseudo/Data/Seqdata/CardosoKaessmann2019NatureRNAseq"
for (Species in c("Human","Mouse")) {
  pseucod = read.csv(paste0("~/Pseudo/Data/Ref/",Species,"/pseudo.coding.closest.bed"),
                     header = F,sep = "\t",stringsAsFactors = F)
  p = list()
  for (Tissue in c("Brain","Cerebellum","Heart","Kidney","Liver","Ovary","Testis")){
    ddpseudo = read.csv(file.path("~/Pseudo/Result",Species,"Savedata/DDG",paste0(Tissue,"_ddg.csv")),
                        header = T,sep = ",",stringsAsFactors = F)
    
    b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",Species,paste0("gene",Species,".bed")), 
                  header = FALSE,sep = "\t",stringsAsFactors = F)
    ddcoding = fread(file.path(wd,Species,paste0(Species,".gene.development.annotation.csv"))) %>% 
      as.data.frame()
    ddcoding %<>% dplyr::select(.,c(ends_with("ID"),ends_with("DDG"))) %>% 
      dplyr::select(.,c(ends_with("ID"),starts_with(Tissue)))
    ddcoding[is.na(ddcoding)]=0
    if (Tissue != "Testis") {
      colnames(ddcoding)[2] = "DDG"
      ddcoding %<>% dplyr::filter(.,DDG >0)
    }else{
      ddcoding = ddcoding[apply(ddcoding[,-1], 1, function(x)sum(x >0) >=1 ),]
      ddcoding$type = b[match(ddcoding[,1],b$V4),5]
      ddcoding %<>% dplyr::filter(.,type==" protein_coding")
    }
    for (i in 1:nrow(pseucod)) {
      if (pseucod[i,4] %in% ddpseudo[,1]){
        pseucod[i,"pseudo"] = 1
      }else{
        pseucod[i,"pseudo"] = 0
      }
    }
    for (i in 1:nrow(pseucod)) {
      if (pseucod[i,11] %in% ddcoding[,1]){
        pseucod[i,"coding"] = 1
      }else{
        pseucod[i,"coding"] = 0
      }
    }
    table(pseucod$pseudo,pseucod$coding) %>% fisher.test(.,alternative = "greater")
    tapply(abs(pseucod$V15), paste0(pseucod$pseudo,pseucod$coding), function(x)median(x,na.rm = T))
    
    pseucod$type = paste0(pseucod$pseudo,pseucod$coding) %>% factor(.,levels = c("11","10","01","00"))
    my_comparisons = list(c("11","10"),c("11","01"),c("11","00"))
    p[[Tissue]] = na.omit(pseucod) %>% ggplot(.,aes(type,log10(abs(V15)),fill=type))+theme_classic()+
      geom_boxplot(notch = T,outlier.colour = "white")+ 
      theme(axis.title.x=element_blank(),#legend.position='none',
            axis.text.x = element_text(size=18),axis.title.y=element_text(size=20),
            axis.text.y = element_text(size=18))+ylab("Distance (log10 bp)")+
      #stat_compare_means(comparisons = my_comparisons)+
      coord_cartesian(ylim = c(2,7))+
      scale_fill_manual(values = c("#b2cee4",rep("#d2d2d2",3)))+
      theme(legend.position='none')+ggtitle(Tissue)+
      theme(plot.title = element_text(hjust = 0.5))
  }
  ggsave(p$Brain+p$Cerebellum+p$Heart+p$Kidney+p$Liver+p$Ovary+p$Testis,
         filename = file.path("~/Pseudo/Result",Species,"Picture",
                              paste0(Species,".Distance.pseu.closest.coding.pdf")),
         device = "pdf",width = 14,height = 9)
}


## Correlation
rm(list = ls())
Species = "Human"
wd = "/home/qians/Pseudo/Data/Seqdata/CardosoKaessmann2019NatureRNAseq"
a <- fread(file.path(wd,Species,paste0(Species,".allgene.tpm.txt"))) %>% as.data.frame()
a[1:3,1:3]
#a <- dplyr::filter(a,!grepl("_PAR_Y",Gene))
#a$V1 <- lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
row.names(a) =a$V1
a =a[,-1]
pseucod = read.csv(paste0("~/Pseudo/Data/Ref/",Species,"/pseudo.coding.closest.bed"),
                   header = F,sep = "\t",stringsAsFactors = F)
express.r = c(); express.p = c()
for (i in 1:nrow(pseucod)) {
  if (pseucod[i,4] %in% row.names(a) & pseucod[i,11] %in% row.names(a)) {
    express.cor = cor.test(as.numeric(a[pseucod[i,4],]),as.numeric(a[pseucod[i,11],]))
    express.r = c(express.r,as.numeric(express.cor$estimate))
    express.p = c(express.p,as.numeric(express.cor$p.value))
  }else{
    express.r = c(express.r,0)
    express.p = c(express.p,1)
  }
}
data = data.frame(Observed=express.r)

c = pseucod[,c(4,11)]
c$disturb = sample(1:nrow(pseucod),size = nrow(pseucod),replace = F)
express.r = c(); express.p = c()
for (i in 1:nrow(c)) {
  if (c[i,1] %in% row.names(a) & c[c[i,3],2] %in% row.names(a)) {
    express.cor = cor.test(as.numeric(a[c[i,1],]),as.numeric(a[c[c[i,3],2],]))
    express.r = c(express.r,as.numeric(express.cor$estimate))
    express.p = c(express.p,as.numeric(express.cor$p.value))
  }else{
    express.r = c(express.r,0)
    express.p = c(express.p,1)
  }
}

data$Simulated = express.r
data[is.na(data)]=0

with(data,wilcox.test(Observed,Simulated,alternative = "greater"))
ggplot(data, aes(x=x) ) +
# Top
geom_density( aes(x = Observed) ) +
geom_label( aes(x=0.6, y=4, label="Observed")) +
geom_vline(xintercept = median(abs(data$Observed)),lty = 2)+
# Bottom
geom_density( aes(x = Simulated)) +
geom_label( aes(x=0.6, y=-4, label="Simulated")) +
geom_vline(xintercept = median(abs(data$Simulated)),lty = 2) +
xlab("Correlation coefficient")+ ylab("Density")+
theme_classic() + theme(
axis.title.x=element_text(size=20),#legend.position='none',
axis.text.x = element_text(size=18),
axis.title.y=element_text(size=20),
axis.text.y = element_text(size=18))
