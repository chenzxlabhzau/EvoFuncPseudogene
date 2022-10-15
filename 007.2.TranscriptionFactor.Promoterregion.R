rm(list = ls())
for (S in c("Human","Mouse")) {
  a = read.csv(file = file.path("~/Pseudo/Data/Ref/",S,paste0(S,".nonreExonlength.2.bed")),header = FALSE,sep = "\t",stringsAsFactors = FALSE)
  a1 = a[a$V5=="+",]
  a2 = a[a$V5=="-",]
  t1 = tapply(a1$V2,a1$V4,min) %>% as.data.frame()
  t2 = tapply(a2$V3,a2$V4,max) %>% as.data.frame()
  t1$tss = t1$. - 2000
  t1$tes = t1$. + 1000
  t2$tss = t2$. - 1000 
  t2$tes = t2$. + 2000
  #Promoter
  promoter = rbind(t1,t2)
  promoter[,c("chr","strand")] = c(a[match(row.names(promoter),a$V4),c(1,5)])
  promoter$name = row.names(promoter)
  promoter$chr = paste0("chr",promoter$chr)
  promoter[promoter$tss<0,"tss"]=0
  TSS = promoter
  promoter = promoter[,c("chr","tss","tes","name","strand")]
  write.table(promoter,file = file.path("~/Pseudo/Data/Ref",S,paste0(S,".promoter.bed")),
              quote = FALSE,row.names = FALSE,sep = "\t",col.names = FALSE)
  #TSS
  TSS$point = 0
  for (i in 1:nrow(TSS)) {
    if (TSS[i,5] == "+") {
      TSS[i,"point"] = TSS[i,"."] + 1
    } else {
      TSS[i,"point"] = TSS[i,"."] - 1
    }
  }
  TSS$point1 = apply(TSS[,c(".","point")],1,min)
  TSS$point2 = apply(TSS[,c(".","point")],1,max)
  TSS = TSS[,c("chr","point1","point2","name","strand")]
  TSS %<>% dplyr::filter(.,chr %in% c(paste0("chr",1:100),"chrX","chrY"))

  write.table(TSS,file = file.path("~/Pseudo/Data/Ref",S,paste0(S,".TSS.bed")),
              quote = FALSE,row.names = FALSE,sep = "\t",col.names = FALSE)
  
  b = read.csv(paste0("~/Pseudo/Data/Ref/",S,"/","gene",S,".bed"),header = FALSE,sep = "\t")
  ddg = fread(paste0("~/Pseudo/Result/",S,"/Savedata/DDG/all.ddg.csv"),header = FALSE) %>% as.data.frame()
  TSS$type = b[match(TSS$name,b$V4),5]
  write.table(TSS[grepl("coding",TSS$type),1:3],file = file.path("~/Pseudo/Data/Ref",S,paste0(S,".TSS.Coding.bed")),
              quote = FALSE,row.names = FALSE,sep = "\t",col.names = FALSE)

  write.table(TSS[grepl("lncRNA",TSS$type),1:3],file = file.path("~/Pseudo/Data/Ref",S,paste0(S,".TSS.lncRNA.bed")),
              quote = FALSE,row.names = FALSE,sep = "\t",col.names = FALSE)
  
  write.table(TSS[grepl("pseu",TSS$type) & TSS$name %in% ddg$V1,1:3],
              file = file.path("~/Pseudo/Data/Ref",S,paste0(S,".TSS.DynamicPseu.bed")),
              quote = FALSE,row.names = FALSE,sep = "\t",col.names = FALSE)
  
  write.table(TSS[grepl("pseu",TSS$type) & !TSS$name %in% ddg$V1,1:3],
              file = file.path("~/Pseudo/Data/Ref",S,paste0(S,".TSS.NodynamicPseu.bed")),
              quote = FALSE,row.names = FALSE,sep = "\t",col.names = FALSE)

  }

