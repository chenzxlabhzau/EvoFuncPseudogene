############
#Human
############
if (TRUE) { # date age with filter paramters
  rm(list = ls())
  a = read.csv("~/Pseudo/Data/Ref/Human/geneHuman.pseudogene.bed",header = F,sep = "\t",stringsAsFactors = F)
  gene.table = matrix(1, nrow = nrow(a),ncol = 13) %>% as.data.frame()
  colnames(gene.table) = c("Human","Chimp","Rhesus","Mouse","Rabbit","Rat","Dog","Opossum","Platypus","Chicken","XTropicalis","Zebrafish","Max")
  gene.table$Human = a$V4
  
  #assign age with score 
  for (n in 2:(ncol(gene.table)-1)) {
    species = colnames(gene.table)[n]
    align = read.csv(paste0("/opt/qians/Pseudogene/Data/UCSC.maf/Human/",species,"/","Human.",species,".Human.pseudogene.exon.alignment.bed"),
                     header = F,sep = "\t",stringsAsFactors = F)
    if (n<=3) {
      gene.table[gene.table$Human %in% align$V4,species]=n
    }else if (n>=4&n<=6) {
      gene.table[gene.table$Human %in% align$V4,species]=4
    }else if (n>=7) {
      gene.table[gene.table$Human %in% align$V4,species]=n-2
    }
  }
  #determine arise time
  gene.table$Max = apply(gene.table[,-c(1,ncol(gene.table))],1,max)
  #Number of support orthologs
  gene.table$SupportNumber = apply(gene.table[,-c(1,ncol(gene.table))],1,function(x)sum(x>1))
  #percent overlap with repeats
  repeats = fread("~/Pseudo/Data/Ref/Human/TE/Human.pseudogene.exon.repeat.txt") %>% 
    separate(.,col = V4,into = c("gene","exon","length"))
  exon.repeat = tapply(as.numeric(repeats$V11), repeats$exon,sum) %>% as.data.frame()
  exon.repeat[,c("gene","length")] = repeats[match(row.names(exon.repeat),repeats$exon),c("gene","length")]
  gene.exon.length = tapply(as.numeric(exon.repeat$length), exon.repeat$gene, sum) %>% as.data.frame()
  gene.repeat.length = tapply(as.numeric(exon.repeat$.), exon.repeat$gene, sum) %>% as.data.frame()
  gene.exon.length$repeatbp = gene.repeat.length[match(row.names(gene.exon.length),row.names(gene.repeat.length)),1]
  gene.exon.length$percentage = gene.exon.length$repeatbp/gene.exon.length$.
  
  gene.table$percentage = gene.exon.length[match(gene.table$Human,row.names(gene.exon.length)),"percentage"]
  gene.table[is.na(gene.table)]=0
  #no Y
  gene.table$chr = a[match(gene.table$Human,a$V4),1]
  gene.table %<>% dplyr::filter(.,percentage <= 0.7 & chr %in% c(1:100,"X"))
  #convert branch to age
  write.csv(gene.table,file = "/home/qians/Pseudo/Result/Human/Savedata/Human.gene.branch.csv",quote = FALSE,row.names = FALSE)
  #details
  with(gene.table[gene.table$Max==3,],table(Chimp,Rhesus))
}

if (TRUE) { #Abundance ~ age
  rm(list = ls())
  library(ggpubr)
  a = read.csv("/home/qians/Pseudo/Result/Human/Savedata/Human.gene.branch.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
  #convert branch to age
  freq = table(a$Max) %>% as.data.frame()
  freq$age = c(-6.65, -22.79, -60.38, -6.64, -62.14, -18.33, -134.98, -39.85, -83.57,-15)
  for (i in rev(2:10)) {
    freq[i,3] = sum(freq[1:(i-1),3]) + freq[i,3]/2
  }
  freq[1,3] = freq[1,3] / 2
  cor = cor.test(as.numeric(freq$age),freq$Freq,method = "spearman")
  p = ggplot(freq,aes(age,log10(Freq)))+geom_point()+geom_smooth(formula = y ~ x, method = "lm",aes(group=" "),se = FALSE)+
    theme_bw()+xlab("Million year")+ylab("Abundance (log10)")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
          axis.text.y = element_text(size=14),axis.title.y = element_text(size=16)#,axis.ticks.y = element_blank(),axis.ticks.x = element_blank()
          )+
    annotate(geom="text", x=-350, y=3.5, size=6,
             label=paste0("rho=",round(cor$estimate,2),",","\n","P=",round(cor$p.value,3)))
  
  p = ggplot(freq,aes(age,log10(Freq)))+geom_point(size=2.5)+geom_smooth(formula = y ~ x, method = "lm",aes(group=" "))+
    geom_line()+
    theme_bw()+xlab("Million year")+ylab("Abundance (log10)")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
          axis.text.y = element_text(size=14),axis.title.y = element_text(size=16)#,axis.ticks.y = element_blank(),axis.ticks.x = element_blank()
    )+
    annotate(geom="text", x=-350, y=3.5, size=6,
             label=paste0("rho=",round(cor$estimate,2),",","\n","P=",round(cor$p.value,3)))
  
  ggsave(p,filename = "/home/qians/Pseudo/Result/Human/Picture/Human.abundance.age.pdf",
         width = 4,height = 4)
  a$age = freq[match(a$Max,freq$Var1),"age"]
  write.csv(a,file = "/home/qians/Pseudo/Result/Human/Savedata/Human.gene.age.csv",quote = FALSE,row.names = FALSE)
  
  ##Only processed 
  S = "Human"
  b = read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
           header = FALSE,sep = "\t",stringsAsFactors = F)
  a$type = b[match(a$Human,b$V4),5]
  a %<>% dplyr::filter(.,grepl("processed",type)) %<>% dplyr::filter(.,!grepl("unprocessed",type))
  freq = table(a$Max) %>% as.data.frame()
  freq$age = c(-6.65, -22.79, -60.38, -6.64, -62.14, -18.33, -134.98, -39.85, -83.57,-15)
  for (i in rev(2:10)) {
    freq[i,3] = sum(freq[1:(i-1),3]) + freq[i,3]/2
  }
  freq[1,3] = freq[1,3] / 2
  cor = cor.test(as.numeric(freq$age),freq$Freq,method = "spearman")
  p = ggplot(freq,aes(age,log10(Freq)))+geom_point(size=2.5)+geom_smooth(formula = y ~ x, method = "lm",aes(group=" "))+
    geom_line()+
    theme_bw()+xlab("Million year")+ylab("Abundance (log10)")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
          axis.text.y = element_text(size=14),axis.title.y = element_text(size=16)#,axis.ticks.y = element_blank(),axis.ticks.x = element_blank()
    )+
    annotate(geom="text", x=-350, y=3.5, size=6,
             label=paste0("rho=",round(cor$estimate,2),",","\n","P=",round(cor$p.value,3)))
  ggsave(p,filename = "/home/qians/Pseudo/Result/Human/Picture/Human.abundance.age.processed.pdf",
         width = 4,height = 4)
  
  if (FALSE) {
    freq$Var1 %<>% rev()
    cor = cor.test(as.numeric(freq$Var1),freq$Freq,method = "spearman")
    ggplot(freq,aes(Var1,log10(Freq)))+geom_point()+geom_smooth(formula = y ~ x, method = "lm",aes(group=" "),se = FALSE)+
      theme_bw()+xlab("From young to old")+#ylab("Abundance (log10)")+
      theme(panel.grid.major =element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text.x = element_blank(),axis.title.x = element_text(size=16),
            axis.text.y = element_blank(),axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank())+
      annotate(geom="text", x=3, y=3.5, size=8,
               label=paste0("rho=",round(cor$estimate,2),",","\n","P=",round(cor$p.value,3)))
    }
  }

############
#Mouse
############
rm(list = ls())
a = read.csv("~/Pseudo/Data/Ref/Mouse/geneMouse.pseudogene.bed",header = F,sep = "\t",stringsAsFactors = F)
gene.table = matrix(1, nrow = nrow(a),ncol = 13) %>% as.data.frame()
colnames(gene.table) = c("Mouse","Rat","Rabbit","Human","Chimp","Rhesus","Dog","Opossum","Platypus","Chicken","XTropicalis","Zebrafish","Max")
gene.table$Mouse = a$V4
#assign age with score 
for (n in 2:(ncol(gene.table)-1)) {
  species = colnames(gene.table)[n]
  align = read.csv(paste0("/opt/qians/Pseudogene/Data/UCSC.maf/Mouse/",species,"/","Mouse.",species,".Mouse.pseudogene.exon.alignment.bed"),
                   header = F,sep = "\t",stringsAsFactors = F)
  if (n<=3) {
    gene.table[gene.table$Mouse %in% align$V4,species]=n
  }else if (n>=4&n<=6) {
    gene.table[gene.table$Mouse %in% align$V4,species]=4
  }else if (n>=7) {
    gene.table[gene.table$Mouse %in% align$V4,species]=n-2
  }
}
#determine arise time
gene.table$Max = apply(gene.table[,-c(1,ncol(gene.table))],1,max)
#Number of support orthologs
gene.table$SupportNumber = apply(gene.table[,-c(1,ncol(gene.table))],1,function(x)sum(x>1))
#percent overlap with repeats
repeats = fread("~/Pseudo/Data/Ref/Mouse/TE/Mouse.pseudogene.exon.repeat.txt") %>% 
  separate(.,col = V4,into = c("gene","exon","length"))
exon.repeat = tapply(as.numeric(repeats$V11), repeats$exon,sum) %>% as.data.frame()
exon.repeat[,c("gene","length")] = repeats[match(row.names(exon.repeat),repeats$exon),c("gene","length")]
gene.exon.length = tapply(as.numeric(exon.repeat$length), exon.repeat$gene, sum) %>% as.data.frame()
gene.repeat.length = tapply(as.numeric(exon.repeat$.), exon.repeat$gene, sum) %>% as.data.frame()
gene.exon.length$repeatbp = gene.repeat.length[match(row.names(gene.exon.length),row.names(gene.repeat.length)),1]
gene.exon.length$percentage = gene.exon.length$repeatbp/gene.exon.length$.

gene.table$percentage = gene.exon.length[match(gene.table$Mouse,row.names(gene.exon.length)),"percentage"]
gene.table[is.na(gene.table)]=0
#no Y
gene.table$chr = a[match(gene.table$Mouse,a$V4),1]
gene.table %<>% dplyr::filter(.,percentage <= 0.7 & chr %in% c(1:100,"X"))
write.csv(gene.table,file = "/home/qians/Pseudo/Result/Mouse/Savedata/Mouse.gene.branch.csv",quote = FALSE,row.names = FALSE)
#details
with(gene.table[gene.table$Max==3,],table(Rat,Rabbit))

if (TRUE) { #Abundance ~ age
  rm(list = ls())
  a = read.csv("/home/qians/Pseudo/Result/Mouse/Savedata/Mouse.gene.branch.csv",header = TRUE,sep = ",",stringsAsFactors = FALSE)
  #convert branch to age
  freq = table(a$Max) %>% as.data.frame()
  freq$age = c(-20.89, -61.25, -7.68, -6.64, -62.14, -18.33, -134.98, -39.85, -83.57,-15)
  for (i in rev(2:10)) {
    freq[i,3] = sum(freq[1:(i-1),3]) + freq[i,3]/2
  }
  freq[1,3] = freq[1,3] / 2
  cor = cor.test(as.numeric(freq$age),freq$Freq,method = "spearman")
  p = ggplot(freq,aes(age,log10(Freq)))+geom_point()+geom_smooth(formula = y ~ x, method = "lm",aes(group=" "),se = FALSE)+
    theme_bw()+xlab("Million year")+ylab("Abundance (log10)")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
          axis.text.y = element_text(size=14),axis.title.y = element_text(size=16)#,axis.ticks.y = element_blank(),axis.ticks.x = element_blank()
    )+
    annotate(geom="text", x=-350, y=3.5, size=6,
             label=paste0("rho=",round(cor$estimate,2),",","\n","P=",round(cor$p.value,3)))
  
  p = ggplot(freq,aes(age,log10(Freq)))+geom_point(size=2.5)+geom_smooth(formula = y ~ x, method = "lm",aes(group=" "))+
    geom_line()+
    theme_bw()+xlab("Million year")+ylab("Abundance (log10)")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
          axis.text.y = element_text(size=14),axis.title.y = element_text(size=16)#,axis.ticks.y = element_blank(),axis.ticks.x = element_blank()
    )+
    annotate(geom="text", x=-350, y=3.5, size=6,
             label=paste0("rho=",round(cor$estimate,2),",","\n","P=",round(cor$p.value,3)))
  
  ggsave(p,filename = "/home/qians/Pseudo/Result/Mouse/Picture/Mouse.abundance.age.pdf",
         width = 4,height = 4)
  a$age = freq[match(a$Max,freq$Var1),"age"]
  write.csv(a,file = "/home/qians/Pseudo/Result/Mouse/Savedata/Mouse.gene.age.csv",quote = FALSE,row.names = FALSE)
  
  ##Only processed 
  S = "Mouse"
  b = read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
               header = FALSE,sep = "\t",stringsAsFactors = F)
  a$type = b[match(a$Mouse,b$V4),5]
  a %<>% dplyr::filter(.,grepl("processed",type)) %<>% dplyr::filter(.,!grepl("unprocessed",type))
  freq = table(a$Max) %>% as.data.frame()
  freq$age = c(-20.89, -61.25, -7.68, -6.64, -62.14, -18.33, -134.98, -39.85, -83.57,-15)
  for (i in rev(2:10)) {
    freq[i,3] = sum(freq[1:(i-1),3]) + freq[i,3]/2
  }
  freq[1,3] = freq[1,3] / 2
  cor = cor.test(as.numeric(freq$age),freq$Freq,method = "spearman")
  p = ggplot(freq,aes(age,log10(Freq)))+geom_point(size=2.5)+geom_smooth(formula = y ~ x, method = "lm",aes(group=" "))+
    geom_line()+
    theme_bw()+xlab("Million year")+ylab("Abundance (log10)")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size=14),axis.title.x = element_text(size=16),
          axis.text.y = element_text(size=14),axis.title.y = element_text(size=16)#,axis.ticks.y = element_blank(),axis.ticks.x = element_blank()
    )+
    annotate(geom="text", x=-350, y=3.5, size=6,
             label=paste0("rho=",round(cor$estimate,2),",","\n","P=",round(cor$p.value,3)))
  ggsave(p,filename = "/home/qians/Pseudo/Result/Mouse/Picture/Mouse.abundance.age.processed.pdf",
         width = 4,height = 4)
}