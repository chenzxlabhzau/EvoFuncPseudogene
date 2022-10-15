rm(list = ls());gc();rm(list = ls())
Num = "015.3."

all = c(
  "ACC",
  "BLCA",
  "BRCA",
  "CESC",
  "CHOL",
  "COAD",
  "DLBC",
  "ESCA",
  "GBM",
  "HNSC",
  "KICH",
  "KIRC",
  "KIRP",
  #"LAML",
  "LGG",
  "LIHC",
  "LUAD",
  "LUSC",
  "MESO",
  "OV",
  "PAAD",
  "PCPG",
  "PRAD",
  "READ",
  "SARC",
  "SKCM",
  "STAD",
  "TGCT",
  "THCA",
  "THYM",
  "UCEC",
  "UCS",
  "UVM"
)

library(ImmuLncRNA)
library(estimate)
library(magrittr)

dir_test <- "~/Pseudo/Data/Seqdata/TCGA/ComRNA4TurPur/"
b = fread("~/Pseudo/Data/Ref/Human/geneHuman.bed") %>% as.data.frame()

for (i in 1:length(all)) {
  focal.cancer = paste0("TCGA-",all[i])
  file_test <- paste0(focal.cancer,".ComRNA.csv")
  #1
  res_test <- ImmuLncRNA::turpur.est(dir_test,file_test)
  print(focal.cancer)
  a = fread(paste0("~/Pseudo/Data/Seqdata/TCGA/MergeCount/",focal.cancer,".mergecount.csv")) %>% as.data.frame()
  row.names(a) = a$V1
  a = a[,-1]
  a %<>% dplyr::select(.,ends_with("TP"))
  colnames(a) %<>% lapply(., function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][1]) %>% unlist()
  a = apply(a, 2, function(x)x/sum(x)*10^6) %>% as.data.frame()
  a$type = b[match(row.names(a),b$V4),"V5"]

  pseudo_exp = a[grepl("pseudo",a$type),-ncol(a)] %>% as.matrix()
  pseudo_exp = pseudo_exp[apply(pseudo_exp, 1, sum) >0,]
  mRNA_exp = a[grepl("protein_coding",a$type),-ncol(a)] %>% as.matrix()
  mRNA_exp = mRNA_exp[apply(mRNA_exp, 1, sum) >0,]
  
  turpur_ori = res_test 
  names(turpur_ori) %<>% lapply(., function(x)strsplit(x,split = "_",fixed = TRUE)[[1]][1]) %>% 
    unlist() %>% gsub(".","-",.,fixed=TRUE)
  
  names(turpur_ori)[1:3];colnames(a)[1:3];colnames(mRNA_exp)[1:3];colnames(pseudo_exp)[1:3]
  #2
  test_res <- par.cor(mRNA_exp,pseudo_exp,turpur_ori,adjusted=FALSE)
  test_res$pcor.value[1:3,1:3]
  
  str(pathways)
  #3
  k=0.995
  test_res <- immu.LncRNA(mRNA_exp,pseudo_exp,adjusted=FALSE,turpur_ori,pathways,k)
  test_res$sig_pairs[1:2,] # showing the immu-lncRNA pairs
  test_res$fgseaRes_all[1:2,] # showing the GSEA results
  table( abs(test_res$fgseaRes_all$sigValue)>k,test_res$fgseaRes_all$padj<0.05)
  write.table(test_res$fgseaRes_all[,-9],file = paste0("/home/qians/Pseudo/Data/Seqdata/TCGA/ImmunoPseudo/",
              focal.cancer,".ImmuPseudo.csv"),quote = FALSE,sep = "\t",row.names = FALSE)
}

