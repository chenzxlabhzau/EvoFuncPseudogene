rm(list = ls());gc();rm(list = ls())
Num = "S003"

library(data.table)
library(magrittr)
library(stringr)
library(dplyr)
library(DESeq2)
  
wd = "/home/qians/Pseudo/Data/Seqdata/Illumina"
Species = "Human"
a <- fread(file.path(wd,Species,"Kaessmann.Mergecount.changed.txt")) %>% as.data.frame()
coldata <- read.csv(file = file.path("~/Pseudo/Result",Species,"Savedata","coldata.csv"),
                    row.names = 1, sep = ",")
coldata$name %<>% gsub("KidneyTestis","Kidney",.)
coldata$time <- lapply(coldata$stage2,function(x){substr(x,1,str_length(x)-3)}) %>% 
  unlist() %>% as.numeric()

if (all(row.names(coldata) %in% colnames(a)[-1])) {
  a <- dplyr::filter(a,!grepl("_PAR_Y",V1))
  row.names(a) = lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
  a = a[,-1]
  a = a[,row.names(coldata)]
}
b = read.csv(paste0("/home/qians/Pseudo/Data/Ref/",Species,"/","gene",Species,".bed"),
             header = F,sep = "\t",stringsAsFactors = FALSE)
b %<>% dplyr::filter(., grepl("pseudogene|protein_coding|lncRNA",V5) & V1 %in% c(1:100,"X"))
b[b$V1!="X","V1"]="A"
b$V5 %<>% gsub(" ","",.)
b[b$V5!="protein_coding"&b$V5!="lncRNA","V5"] ="Pseudogene"
b[b$V5=="protein_coding","V5"] ="Coding"

coldata$condition %<>% gsub("Testis","Gonad",.) %>% gsub("Ovary","Gonad",.)
c = c()
for ( Tissue in unique(coldata$condition) ) {
  coldata2 = dplyr::filter(coldata, condition == Tissue )
  #coldata2$stage4 %<>% factor(.,levels = deve[Tissue,][deve[Tissue,]!=" "] )
  onetissue = a[,row.names(coldata2)] #one tissue matrix
  if (all(row.names(coldata2) == colnames(onetissue))) {
    dds = DESeq2::DESeqDataSetFromMatrix(onetissue,
                                         colData = coldata2,
                                         design = ~sex)
  }
  dds <- DESeq(dds)
  resultsNames(dds) # lists the coefficients
  
  res = results(dds,contrast = c("sex","Male","Female")) %>% as.data.frame() %>%
    #dplyr::filter(.,abs(log2FoldChange) >=0.5 & padj <= 0.05) %>% 
    dplyr::select(.,c(log2FoldChange,padj))
  res[,c("type","chr")] = b[match(row.names(res),b$V4),c(5,1)]
  res %<>% na.omit #%>% dplyr::filter(., grepl("pseudogene|protein_coding",type) & chr %in% c(1:100,"X")) 
  res$bias = "Unbiased"
  res[res$log2FoldChange>=0.5 & res$padj <= 0.05, "bias"] = "Male-biased"
  res[res$log2FoldChange<= -0.5 & res$padj <= 0.05, "bias"] = "Female-biased"
  
  df = as.data.frame(t(apply(as.matrix(table(res$bias,res$type)), 1, function(x)x/sum(x))))
  df$bias = row.names(df)
  df$tissue = Tissue
  df %<>% pivot_longer(.,cols=1:3)
  df$name %<>% gsub("lncRNA","LncRNA",.) %>% factor(.,levels = rev(unique(.))) 
  c = rbind(c,df)
}
c$tissue %<>% factor(.,levels = c("Brain","Cerebellum","Heart","Liver","Kidney","Gonad"))
  p = ggplot(c, aes(x = 1, y = value, fill = name)) +
    geom_bar(stat = "identity", width = 1) +    
    coord_polar(theta = "y")+facet_grid(bias~tissue,switch = "y")+theme_classic()+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
          axis.title.y = element_blank(),axis.text.y = element_blank(),
          axis.line.x = element_blank(), axis.line.y = element_blank(),
          axis.ticks = element_blank(),
          strip.text.x = element_text(size = 12),strip.text.y = element_text(size = 12,angle=0),
          strip.background = element_blank(),strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom",legend.direction = "horizontal",legend.title = element_blank())+
    geom_text(aes(label = round(value*100,2)), position = position_stack(vjust = 0.5))+
    scale_fill_manual(values =c("#fc9272","#9ecae1","#d95f0d"))
  ggsave(p, filename = paste0("/home/qians/Pseudo/Result/",Species,"/Picture/",Num,".Sexbias.tissue.pdf"),
         width = 9,height = 6)
    
  species = Species
  age = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),header = TRUE,sep = ",",stringsAsFactors = FALSE)
  age$bias = res[match(age$Human,row.names(res)),"bias"]
  age[is.na(age$bias),"bias"]="Unbiased"
  cor.test(as.numeric(table(age$bias=="Male-biased",age$age)[2,]/table(age$age)),1:10,method = "spearman")
  with(age[age$chr!="X",],table(bias=="Male-biased",age))[2,]/with(age[age$chr!="X",],table(age))
  
  ratiomale = as.data.frame(table(age$bias,age$age)[2,]/table(age$age),stringsAsFactors=FALSE)
  ratiomale$bias = "Male-biased"
  ratiomale$Var1 %<>% as.numeric()
  cor = cor.test(ratiomale$Var1,ratiomale$Freq,method = "spearman")
  p = ggplot(ratiomale,aes(Var1,Freq))+geom_point()+geom_smooth(method = "lm")+theme_bw()+
    xlab("Million year")+ylab("Proportion of male-biased pseudogene (%)")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size=10),axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=10),axis.title.y = element_text(size=12))+
    annotate(geom="text", x=-400, y=0.05, size=6,
             label=paste0("rho=",round(cor$estimate,2),",","\n","P=",round(cor$p.value,3)))
  ggsave(p, filename = paste0("/home/qians/Pseudo/Result/Human/Picture/",Num,".Human.Proportion.Malebias.Gonad.Age.pdf"),
         width = 5,height = 4)
  
  ratiounbias = as.data.frame(table(age$bias,age$age)[3,]/table(age$age),stringsAsFactors=FALSE)
  ratiounbias$bias = "Unbiased"
  ratiounbias$Var1 %<>% as.numeric()
  cor = cor.test(ratiounbias$Var1,ratiounbias$Freq,method = "spearman")
  p = ggplot(ratiounbias,aes(Var1,Freq))+geom_point()+geom_smooth(method = "lm")+theme_bw()+
    xlab("Million year")+ylab("Proportion of unbiased pseudogene (%)")+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size=10),axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=10),axis.title.y = element_text(size=12))+
    annotate(geom="text", x=-400, y=0.9, size=6,
             label=paste0("rho=",round(cor$estimate,2),",","\n","P=",round(cor$p.value,4)))
  ggsave(p, filename = paste0("/home/qians/Pseudo/Result/Human/Picture/",Num,".Human.Proportion.Unbias.Gonad.Age.pdf"),
         width = 5,height = 4)
  
  
  if (FALSE) {
    species="Human"
    age = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),header = TRUE,sep = ",",stringsAsFactors = FALSE)
    age$bias = res[match(age$Human,row.names(res)),"bias"]
    cor.test(as.numeric(table(age$bias=="Female-biased",age$age)[2,]/table(age$age)),1:10,method = "spearman")
    with(age[age$chr=="X",],table(bias=="Male-biased",age))[2,]/with(age[age$chr=="X",],table(age))
    
    table(res$bias,paste(res$type,res$chr))
    table(res$bias,paste(res$type,res$chr))[c(1,2),3:4]
    
    
    sexbias = table(paste(res$type,res$chr,sep = "_"),res$log2FoldChange>0) %>% as.data.frame()
    sexbias$Var2 %<>% gsub("FALSE","Female-biased",.) %>% gsub("TRUE","Male-biased",.) 
    all = table(paste(b$V5,b$V1,sep="_")) %>% as.data.frame()
    sexbias$all = all[match(sexbias$Var1,all$Var1),2]
    sexbias$not = sexbias$all - sexbias$Freq
    c = data.frame(Type=unique(paste(sexbias$Var1,sexbias$Var2,sep = "_") %>% gsub("_A","",.) %>% gsub("_X","",.)))
    for (i in 1:4) {
      p = fisher.test(sexbias[c(2*i-1,2*i),c(3,5)])$p.value
      odd = fisher.test(sexbias[c(2*i-1,2*i),c(3,5)])$estimate %>% as.numeric()
      c[i,2] = p 
      c[i,3] = odd 
    }
    c$tissue = Tissue
    df = rbind(df,c)
  }
  df$type = lapply(df$Type,function(x)strsplit(x,split = "_")[[1]][1]) %>% unlist()
  df$sex = lapply(df$Type,function(x)strsplit(x,split = "_")[[1]][2]) %>% unlist()
  ggplot(df,aes())
  }

