rm(list = ls());gc();rm(list = ls())
Num = "015.5."

#infor downloaded from TIMER2: http://timer.cistrome.org/
#wget http://timer.cistrome.org/infiltration_estimation_for_tcga.csv.gz

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

infil = fread("/home/qians/Pseudo/Data/Seqdata/TCGA/Infiltration/infiltration_estimation_for_tcga.csv")%>% as.data.frame()
row.names(infil) = infil$cell_type
infil = infil[2:7]
b = fread("~/Pseudo/Data/Ref/Human/geneHuman.bed") %>% as.data.frame()

for (i in 1:length(all)) {
  focal.cancer = paste0("TCGA-",all[i]) 
  a = fread(file = paste0("~/Pseudo/Data/Seqdata/TCGA/MergeCount/",focal.cancer,".mergecount.csv")) %>% as.data.frame()
  row.names(a) = a$V1
  a = a[,-1]
  a %<>% dplyr::select(.,ends_with("TP"))
  colnames(a) %<>% substr(.,1,15)
  a = apply(a, 2, function(x)x/sum(x)*10^6) %>% as.data.frame()
  a = a[apply(a, 1, sum)>0,]
  a$type = b[match(row.names(a),b$V4),"V5"]
  a = a[grepl("pseudo",a$type),-ncol(a)]
  subinfil = infil[colnames(a),]
  
  df = matrix(0,nrow = nrow(a),ncol = 12) %>% as.data.frame()
  row.names(df)  = row.names(a)
  colnames(df) =  paste(rep(c("B","CD4","CD8","Neutrophil","Macrophage","Dendritic"),each=2),c("R","P"),sep = ".")
  
  df[,seq(1,11,2)] = lapply(1:6, function(j){
        apply(a, 1, function(x)cor.test(x,subinfil[,j],method = "spearman")$estimate) %>% as.numeric()
      })
  df[,seq(2,12,2)] = lapply(1:6, function(j){
    apply(a, 1, function(x)cor.test(x,subinfil[,j],method = "spearman")$p.value) %>% as.numeric()
  })
  
  #if (colnames(a) == row.names(subinfil)) {
  #  for (j in 1:6){
  #    df[,2*j-1] = apply(a, 1, function(x)cor.test(x,subinfil[,j],method = "spearman")$estimate) %>% as.numeric()
  #    df[,2*j] = apply(a, 1, function(x)cor.test(x,subinfil[,j],method = "spearman")$p.value) %>% as.numeric()
  #  }
  #}
  write.table(df,file = paste0("~/Pseudo/Data/Seqdata/TCGA/Infiltration/",focal.cancer,".infiltration.cor.csv"),quote = FALSE,sep = "\t")
}
