rm(list = ls())
library(maSigPro)
S="Mouse"
wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
a <- fread(file.path(wd,S,paste0(S,".pseudogene.tpm.txt"))) %>% as.data.frame()
a[1:3,1:3]
# a <- dplyr::filter(a,!grepl("_PAR_Y",V1))
# a$V1 <- lapply(a$V1,function(x)strsplit(x,".",fixed=T)[[1]][1]) %>% unlist()
# b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
#               header = FALSE,sep = "\t",stringsAsFactors = F)
# a$type <- b[match(a$V1,b$V4),5]
# a <- dplyr::filter(a,grepl("pseudogene",type,ignore.case = T))
row.names(a) <- a$V1
a <- a[,-1]
a %<>% dplyr::select(contains("rep"))
# a <- apply(a, 2, function(x){10^6 * x/sum(x)})
express = fread(file.path("~/Pseudo/Result",S,"Savedata","gene.expressnum.0.3allfpkm.csv")) 
express %<>% dplyr::filter(.,num >= 1)
a <- a[express$gene,]

coldata <- read.csv(file = file.path("~/Pseudo/Result",S,"Savedata","coldata.csv"),
                    row.names = 1, sep = ",")
for (i in 1:nrow(coldata)) {
  if (substr(coldata[i,2],1,1)=="e") {
    coldata[i,6] = substr(coldata[i,2],2,5) %>% as.numeric()
  }
  if (substr(coldata[i,2],2,4)== "dpb") {
    coldata[i,6] = substr(coldata[i,2],1,1) %>% as.numeric() +20 #20 day
  }
  if (substr(coldata[i,2],2,4)== "wpb") {
    coldata[i,6] = substr(coldata[i,2],1,1) %>% as.numeric() * 7 +20 #20 day
  }
}
coldata$V6 <- log2(coldata$V6)
coldata <- coldata[order(coldata$V6),]
for (j in as.character(unique(coldata$condition))) {
  c1 <- dplyr::filter(coldata,condition == j)
  d <- table(c1$V6) %>% as.data.frame()
  c1$Replicates <- rep(1:nrow(d),times=d$Freq)
  if (length(unique(paste0(c1$V6,"_",c1$Replicates))) == length(unique(c1$V6))) {
    row.names(c1) <- c1$name
    c1 <- c1[,c("V6","Replicates")]
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


### Revise 1.
rm(list = ls());gc();rm(list = ls())
for (Species in c("Human","Mouse")) {
  dir = file.path("/home/qians/Pseudo/Result",Species,"Savedata/DDG")
  Files = grep("_ddg.csv$",list.files(dir),value=TRUE)
  filePath <- sapply(Files,function(x){paste(dir,x,sep='/')})
  data <- lapply(filePath, function(x){ fread(x, header=T,sep = ",")})
  a = data[[1]] %>% as.data.frame()
  a$tissue =  gsub("_ddg.csv","",names(data))[1]
  
  for (i in 2:length(data)) {
    b <- data[[i]] %>% as.data.frame()
    b$tissue =  gsub("_ddg.csv","",names(data))[i]
    a <- rbind(a,b)
  }
  fwrite(a,file = paste0("/home/qians/Pseudo/Result/",Species,"/Savedata/DDG/all.ddg.",Species,".csv"))
}


##### compare Pacbio detected and DDPs
rm(list = ls());gc();rm(list = ls())
Species="Mouse"

ddp = fread(paste0("/home/qians/Pseudo/Result/",Species,"/Savedata/DDG/all.ddg.",Species,".csv"))
pacbio = fread(paste0("~/Pseudo/Data/Seqdata/PacBio/mouse_pseudogene.IsoseqDetected.bed"),header = FALSE) %>% as.data.frame()
table(pacbio$V3 %in% unique(ddp$V1))