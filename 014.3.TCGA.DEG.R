rm(list = ls());gc();rm(list = ls())
Num = "014.3."

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

b = read.csv("~/Pseudo/Data/Ref/Human/geneHuman.bed",header = FALSE,sep = "\t")
ddg = fread("~/Pseudo/Result/Human/Savedata/DDG/all.ddg.csv",header = FALSE) %>% as.data.frame()

for (i in 1:length(all)) {
  focal.cancer = paste0("TCGA-",all[i])
  a = fread(paste0("~/Pseudo/Data/Seqdata/TCGA/MergeCount/",focal.cancer,".mergecount.csv")) %>% as.data.frame()
  row.names(a) = a$V1
  a = a[,-1]
  if ( sum(lapply(colnames(a), function(x)strsplit(x,split="-",fixed=TRUE)[[1]][9]) %>%unlist() == "TP")>=5 &
       sum(lapply(colnames(a), function(x)strsplit(x,split="-",fixed=TRUE)[[1]][9]) %>%unlist() == "NT")>=5 ) {
    info = colnames(a) %>% as.data.frame()
    info$type = lapply(info$., function(x)strsplit(x,split="-",fixed=TRUE)[[1]][9]) %>%unlist()
    dds = DESeqDataSetFromMatrix(countData = a,colData = info,design = ~type)
    dds = DESeq(dds)
    resultsNames(dds)
    res = results(dds, contrast = c("type","TP","NT")) %>% as.data.frame()
    res = res[!is.na(res$padj),]
    res = res[abs(res$log2FoldChange)>2&res$padj<0.01,]
    res$type = b[match(row.names(res),b$V4),"V5"]
    res[is.na(res)] = "NA"
    res$dynamic = "NO"
    res[row.names(res) %in% ddg$V1,"dynamic"] = "DDG"
    write.table(res,file = paste0("~/Pseudo/Data/Seqdata/TCGA/DEG/",focal.cancer,".DEG.csv"),quote = FALSE,sep = "\t")
  }else {
    print(focal.cancer)
  }
}

#Download clinic info:
#clin <- GDCquery_clinic("TCGA-ACC", type = "clinical", save.csv = TRUE) 
#Subtype:
#https://www.facingourrisk.org/info/risk-management-and-treatment/cancer-treatment/by-cancer-type/breast/stages-and-subtypes