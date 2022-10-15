rm(list = ls());gc();rm(list = ls())
Num = "014.1."

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
  "LAML",
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

library(TCGAbiolinks)
for (i in 1:length(all)) {
  focal.cancer = paste0("TCGA-",all[i])
  query <- GDCquery(
    project = focal.cancer,
    data.category = "Transcriptome Profiling", 
    data.type = "Gene Expression Quantification",
    workflow.type = "HTSeq - Counts")
  #GDCdownload(query,directory = "/home/qians/Pseudo/Data/Seqdata/TCGA/RawCount/")
  data <- GDCprepare(query,directory = "/home/qians/Pseudo/Data/Seqdata/TCGA/RawCount/")

  # TP,NT,TM,TR,... http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html#Filter_functions
  #meta = colData(data) %>% as.data.frame()
  if ( all(c("TP","NT") %in% unique(data$shortLetterCode)) ) {
    TP = TCGAquery_SampleTypes(data$barcode,"TP")
    NT = TCGAquery_SampleTypes(data$barcode,"NT")
    TP = assay(data)[,TP]
    NT = assay(data)[,NT,drop=FALSE]
    colnames(TP) = paste0(colnames(TP),"_",focal.cancer,"-","TP")
    colnames(NT) = paste0(colnames(NT),"_",focal.cancer,"-","NT")
    mtx = cbind(TP,NT)
    write.table(mtx,file = paste0("~/Pseudo/Data/Seqdata/TCGA/MergeCount/",focal.cancer,".mergecount.csv"),
                quote = FALSE,sep = "\t")
  }
}
#except "LAML",
for (i in c("ACC","DLBC","ESCA","LGG","MESO","OV","TGCT","UCS","UVM")) {
  focal.cancer = paste0("TCGA-",i)
  query <- GDCquery(
    project = focal.cancer,
    data.category = "Transcriptome Profiling", 
    data.type = "Gene Expression Quantification",
    workflow.type = "HTSeq - Counts")
  data <- GDCprepare(query,directory = "/home/qians/Pseudo/Data/Seqdata/TCGA/RawCount/")
  TP = TCGAquery_SampleTypes(data$barcode,"TP")
  TP = assay(data)[,TP]
  colnames(TP) = paste0(colnames(TP),"_",focal.cancer,"-","TP")
  write.table(TP,file = paste0("~/Pseudo/Data/Seqdata/TCGA/MergeCount/",focal.cancer,".mergecount.csv"),
              quote = FALSE,sep = "\t")
}

