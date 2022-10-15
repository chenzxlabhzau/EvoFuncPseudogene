library(Hmisc)
rm(list = ls())
gc()
wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
for (S in c("Human","Mouse")) {
  exp = fread(file.path(wd,S,paste0(S,".allgene.fpkm.txt"))) %>% as.data.frame()
  row.names(exp) = exp$V1
  exp =exp[,-1]
  exp = exp[apply(exp, 1, function(x)sum(x > 0.1) >= 1),]
  exp %<>% dplyr::select(.,!starts_with("Testis")) # exclude testis sample
  b <- read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
                header = FALSE,sep = "\t",stringsAsFactors = F)
  exp$type <- b[match(row.names(exp),b$V4),5]
  
  coding = dplyr::filter(exp,grepl("protein_coding",type,ignore.case = T)) %>% 
    dplyr::select(.,-type)
  coding = coding[apply(coding, 1, function(x)sum(x > 1) >=5),]
  pseu = dplyr::filter(exp,grepl("pseudogene",type,ignore.case = T)) %>% 
    dplyr::select(.,-type)
  pseu = pseu[apply(pseu, 1, function(x)sum(x > 0.1) >=5),]
  
  # if (FALSE) {
  #   compute <- function(x, y, ...) {
  #     t <- lapply(X = x, FUN = function(x.1, y.1 = y, ...) {
  #       correlation <- lapply(X = y.1, cor, x.1, ...)
  #       p.value <- lapply(X = y.1, function(x.2, y.2 = x.1, ...) {
  #         p <- cor.test(x.2, y.2, ...)
  #         p[["p.value"]]
  #       }, ...)
  #       correlation <- do.call(data.frame, correlation)
  #       correlation <- do.call(rbind, correlation)
  #       p.value <- do.call(data.frame, p.value)
  #       p.value <- do.call(rbind, p.value)
  #       result <- data.frame(correlation = correlation[, 1],
  #                            p.value = p.value[, 1])
  #     }, ...)
  #     return(do.call(rbind, t))
  #   }
  # }
  
  # c = c()
  # for (i in 1:2) {
  #   apply(pseu[1:10,], 1, function(x){
  #     cor = cor.test(x,as.numeric(coding[i,]))
  #     r = as.numeric(cor$estimate); p = as.numeric(cor$p.value)
  #     return(c(r,p))
  #   })
  # }
  
  merge_dat <- cbind(t(coding), t(pseu))
  cc<- rcorr(merge_dat, type="pearson")
  trans <- function(object){
    pvalue = object[, !(colnames(object) %in% rownames(pseu))]
    pvalue = pvalue[rownames(object) %in% rownames(pseu),]
    result = data.frame(pseu = rep(rownames(pvalue), ncol(pvalue)),
                        mRNA = rep(colnames(pvalue), each = nrow(pvalue)),
                        values = matrix(pvalue, ncol = 1)[,1])
    return(result)
  }
  PVALUE <- trans(cc$P)
  Correlation <- trans(cc$r)
  PVALUE$ID<- paste(PVALUE$pseu,PVALUE$mRNA,sep = "_")
  Correlation$ID<- paste(Correlation$pseu,Correlation$mRNA,sep = "_")
  if (all(Correlation$ID == PVALUE$ID)) {
    Pair<- cbind(Correlation,PVALUE)
  }
  Pair = Pair[,c(1,2,3,7)]
  colnames(Pair) = c("psue","coding","R","Pvalue")
  fwrite(Pair,file = file.path("~/Pseudo/Result",S,"Savedata","expressionPair.csv"),nThread = 5)
  #write.csv(Pair,file = file.path("~/Pseudo/Result",S,"Savedata","expressionPair.csv"),
  #         quote = FALSE,row.names = FALSE)
}

