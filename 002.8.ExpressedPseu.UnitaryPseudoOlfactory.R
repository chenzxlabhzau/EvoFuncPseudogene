rm(list = ls());gc();rm(list = ls())
Num = "002.8."

exp = fread("~/Pseudo/Result/Human/Savedata/gene.expressnum.1allfpkm.csv") %>% as.data.frame()
S="Human"
b = read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
             header = FALSE,sep = "\t",stringsAsFactors = F)
exp$type = b[match(exp$gene,b$V4),5]

tmp = dplyr::filter(exp,grepl("unitary",type))

id = fread("~/Pseudo/Data/Ref/Human/gene.ENSG.name.2.txt",header = FALSE) %>% as.data.frame()
tmp$id = id[match(tmp$gene,id$V1),2]

vn = read.csv("~/Pseudo/Data/Ref/Human/Unitary.Pseudo.ENSEMBLannota.txt",header = TRUE,sep = "\t")
vn %<>% dplyr::filter(.,grepl("vomeronasal|olfactory",Gene.description))
table(tmp$gene %in% vn$Gene.stable.ID, tmp$num>=1) %>% fisher.test()
