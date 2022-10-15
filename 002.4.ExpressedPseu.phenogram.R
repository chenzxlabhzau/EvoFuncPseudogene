rm(list = ls())
wd = "/home/qians/Pseudo/Data/Seqdata/CardosoKaessmann2019NatureRNAseq"
num = c() 
for (S in c("Human","Mouse")) {
  rpkm = fread(file.path(wd,S,paste0(S,".pseudogene.fpkm.txt"))) %>% as.data.frame()
  row.names(rpkm) = rpkm$V1
  rpkm = rpkm[,-1]
  name = colnames(rpkm) %>% gsub("KidneyTestis","Kidney",.) %>% 
    gsub("Forebrain","Brain",.) %>% 
    gsub("Hindbrain","Cerebellum",.) %>% as.data.frame() 
  #Tissue
  name$Tissue = lapply(name$.,function(x)strsplit(as.character(x),"_",fixed=T)[[1]][1]) %>% unlist()
  rpkm2 = t(apply(rpkm, 1,  function(x) tapply(x, name$Tissue, max))) %>% as.data.frame()
  #rpkm3 = t(apply(rpkm, 1,  function(x) tapply(x, name$Tissue, 
  #                                             function(x)colnames(rpkm)[which.max(x)])))%>%as.data.frame
  rpkm2$max = apply(rpkm2,1,max)
  rpkm2 = rpkm2[order(rpkm2$max,decreasing = T),]
  
  N = 1000
  rpkm2 = rpkm2[1:N,-ncol(rpkm2)]
  
