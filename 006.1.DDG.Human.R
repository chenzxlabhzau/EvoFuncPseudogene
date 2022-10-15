rm(list = ls());gc();rm(list = ls())
library(maSigPro)
S="Human"
wd = "/home/qians/Pseudo/Data/Seqdata/Illumina"
a <- fread(file.path(wd,S,paste0(S,".pseudogene.tpm.txt"))) %>% as.data.frame()
a[1:3,1:3]
#a <- dplyr::filter(a,!grepl("_PAR_Y",V1))
#a$V1 <- lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
#b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
 #             header = FALSE,sep = "\t",stringsAsFactors = F)
#a$type <- b[match(a$V1,b$V4),5]
#a <- dplyr::filter(a,grepl("pseudogene",type,ignore.case = T))
row.names(a) <- a$V1
a <- a[,-1]
#a <- a[,-c(1,ncol(a))]
#a <- apply(a, 2, function(x){10^6 * x/sum(x)})
express = fread(file.path("~/Pseudo/Result",S,"Savedata","gene.expressnum.0.3allfpkm.csv")) 
express %<>% dplyr::filter(.,num >= 1)
a <- a[express$gene,]

coldata <- read.csv(file = file.path("~/Pseudo/Result",S,"Savedata","coldata.csv"),
                    row.names = 1, sep = ",")
coldata$time <- lapply(coldata$stage2,function(x){substr(x,1,str_length(x)-3)}) %>% unlist() %>% as.numeric()
for (i in 1:nrow(coldata)) {
  if (coldata[i,4]=="week") {
    coldata[i,9] = coldata[i,8] * 7
  }
  if (coldata[i,4]=="day") {
    coldata[i,9] = coldata[i,8] * 1 + 280
  }
  if (coldata[i,4]=="month") {
    coldata[i,9] = coldata[i,8] * 30 + 280
  }
  if (coldata[i,4]=="year") {
    coldata[i,9] = coldata[i,8] * 365 + 280
  }
}
coldata$V9 <- log2(coldata$V9)
coldata <- coldata[order(coldata$V9),]
for (j in as.character(unique(coldata$condition))) {
  c1 <- dplyr::filter(coldata,condition == j)
  d <- table(c1$V9) %>% as.data.frame()
  c1$Replicates <- rep(1:nrow(d),times=d$Freq)
  if (length(unique(paste0(c1$V9,"_",c1$Replicates))) == length(unique(c1$V9))) {
    row.names(c1) <- c1$name
    c1 <- c1[,c("V9","Replicates")]
    colnames(c1) <- c("Time","Replicates")
    c1$Control <- 1
    a1 <- a[,row.names(c1)]
    #q1 <- make.design.matrix(c1,degree = 3)
    #fit <- p.vector(a1, q1, Q = 0.05, min.obs = 5)
    #tstep <- T.fit(fit, step.method ="backward", alfa = 0.05)
    #sol <- tstep$sol
    #ddg <- sol[sol$`R-squared` > 0.3,]
    #sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")
    #see.genes(data = sigs$sig.genes$Control,dis = q1$dis, show.fit = T, cluster.method="hclust" ,cluster.data = 1, k = 9,newX11 = F)
    tc.test <- maSigPro(data = a1, edesign = c1, degree = 3,Q = 0.05, min.obs = 5,step.method = "backward", rsq = 0.3)
    d2 <-tc.test$sig.genes$Control$sig.pvalues
    write.csv(d2,file = file.path("/home/qians/Pseudo/Result",S,"Savedata/DDG",paste0(j,"_ddg.csv")),quote = F)
  }
}

# Revise 1 
for (S in c("Human","Mouse")) {
  a = fread(file.path("/home/qians/Pseudo/Result",S,"Savedata/DDG/all.ddg.csv"),header = FALSE) %>% as.data.frame()
  a = dplyr::distinct(a,V1) 
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                               header = FALSE,sep = "\t",stringsAsFactors = F)
  a$type = b[match(a$V1,b$V4),5]
  a$type %<>% gsub("transcribed_","",.) %>% gsub("translated_","",.) %>%
    gsub(" ","",.) %>%
    gsub("unitary_pseudogene","Unitary",.) %>% 
    gsub("polymorphic_pseudogene","Polymorphic",.) %>% 
    gsub("processed_pseudogene","Processed",.)
  a[!a$type %in% c("Unitary","Polymorphic","Processed"),"type"]="Unprocessed"
  print(c(nrow(a),S))
  print(table(a$type))
}
