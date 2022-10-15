rm(list = ls());gc();rm(list = ls())
Num = "015.2."

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
#awk -F "\t" '{print $9"\t"$10}' gene.ENSG.name.txt | sed 's/gene_id //g' | sed 's/gene_name //g' | sed 's/"//g' > gene.ENSG.name.2.txt
id = fread("~/Pseudo/Data/Ref/Human/gene.ENSG.name.2.txt",header = FALSE) %>% as.data.frame()
comgene = fread("~/R/x86_64-pc-linux-gnu-library/4.0/ImmuLncRNA/extdata/mRNA_test.txt") %>% as.data.frame()
#pseudo = fread("~/Pseudo/Data/Ref/Human/geneHuman.pseudogene.bed") %>% as.data.frame()

for (i in 1:length(all)) {
  focal.cancer = paste0("TCGA-",all[i])
  a = fread(paste0("~/Pseudo/Data/Seqdata/TCGA/MergeCount/",focal.cancer,".mergecount.csv")) %>% as.data.frame()
  row.names(a) = a$V1
  a = a[,-1]
  a %<>% dplyr::select(.,ends_with("TP"))
  a = apply(a, 2, function(x)x/sum(x)*10^6) %>% as.data.frame()
  a$symbol = id[match(row.names(a),id$V1),2]
  comgene.express = a[a$symbol %in% comgene$V1,]
  comgene.express = apply(comgene.express[,-ncol(comgene.express)],2, function(x)tapply(x,comgene.express$symbol,mean))
  write.table(comgene.express,file = paste0("~/Pseudo/Data/Seqdata/TCGA/ComRNA4TurPur/",focal.cancer,".ComRNA.csv"),quote = FALSE,sep = "\t")
}
