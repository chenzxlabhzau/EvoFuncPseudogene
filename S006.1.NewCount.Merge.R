# Prepare files -----------------------------------------------------------
#cd ~/ddg_lncRNA/data/nature_ddg/human && mkdir -p ~/ddg_lncRNA/result/human &&  STAR --runThreadN 10 --genomeDir ~/ddg_lncRNA/reference/human  --readFilesCommand zcat --readFilesIn ERR2598067.fastq.gz --outSAMstrandField intronMotif --outFileNamePrefix ~/ddg_lncRNA/result/human/ERR2598067 && samtools view -Sb -q 30 -h -@ 10 ~/ddg_lncRNA/result/human/ERR2598067Aligned.out.sam |samtools sort -@ 10 >~/ddg_lncRNA/result/human/ERR2598067.bam && featureCounts -T 10 -s 2 -t exon -g gene_id -a ~/ddg_lncRNA/reference/human_count.gtf ~/ddg_lncRNA/result/human/ERR2598067.bam -o ~/ddg_lncRNA/result/human/ERR2598067.count && rm ERR2598067.fastq.gz && rm ~/ddg_lncRNA/result/human/ERR2598067Aligned.out.sam
#scp lchen@211.69.141.147:/home/lchen/ddg_lncRNA/result/*.count ./

#cd ~/Pseudo/Data/Seqdata/Illumina/Human
#wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6814/E-MTAB-6814.sdrf.txt
#cd ~/Pseudo/Data/Seqdata/Illumina/Mouse
#wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6798/E-MTAB-6798.sdrf.txt

#cd /NAS/qians/NewCount/Mouse/Ourdata
#

# Merge Kaessmann ---------------------------------------------------------
#* Merge count ------------------------------------
if (FALSE) {# has finished
  rm(list = ls());gc();rm(list = ls())
  Num = "S006.1."
  for (Species in c("Human","Mouse")) {
    dir = paste0("/NAS/qians/NewCount/",Species,"/Kaessmann/")
    Files = grep("count$",list.files(dir),value=TRUE)
    filePath <- sapply(Files,function(x){paste(dir,x,sep='/')})
    data <- lapply(filePath, function(x){ fread(x, header=T,sep = "\t")})
    a = data[[1]] %>% as.data.frame()
    a = a[,c(1,ncol(a))]
    colnames(a)[2] %<>% gsub(".bam","",.)
    
    for (i in 2:length(data)) {
      b <- data[[i]] %>% as.data.frame()
      b = b[,c(1,ncol(b))]
      colnames(b)[2] %<>% gsub(".bam","",.)
      a <- merge(a,b,by="Geneid")
    }
    row.names(a) = a$Geneid
    a = a[,-1]
    colnames(a) %<>% gsub("/home/lchen/ddg_lncRNA/result/human/","",.)
    fwrite(a,file = paste0("~/Pseudo/Data/Seqdata/Illumina/",Species,"/Kaessmann.Mergecount.raw.txt"),
           col.names = TRUE,row.names = TRUE, quote = FALSE,sep = "\t")
    rm(list = ls());gc();rm(list = ls())
  }
}

#* Change count matrix colname ------------------------------------
#** Human ------------------------------------
rm(list = ls());gc();rm(list = ls())
Num = "S006.1."
for (Species in "Human") {
  dir = paste0("~/Pseudo/Data/Seqdata/Illumina/",Species)
  a = grep("sdrf.txt$",list.files(dir),value=TRUE) %>% sapply(.,function(x){paste(dir,x,sep='/')}) %>%
    lapply(., function(x){ fread(x, header=T,sep = "\t")})
  a = a[[1]][,c(30,34)] %>% as.data.frame()
  a$V3 <- lapply(a$`Derived Array Data File`,function(x)strsplit(x,split = paste0(Species,"."))[[1]][2]) %>% unlist()
  a$V4 <- lapply(a$V3,function(x)strsplit(x,split = ".sorted")[[1]][1]) %>% unlist()
  a <- a[order(a$V4),]
  f <- table(a$V4) %>% as.data.frame() 
  f$Var1 <- as.character(f$Var1)
  f <- f[order(f$Var1),]
  q <- c()
  for (i in 1:nrow(f)) {
    q <- c(q,paste0(f[i,1],"_rep",1:(f[i,2])))
  }
  a$name <- q %>% gsub(".","_",.,fixed = TRUE)
  # a$name2 <- lapply(a$name,function(x)strsplit(x,split = "_rep")[[1]][1]) %>% unlist()
  # table(a$V4==a$name2) %>% as.numeric() == length(unique(sort(a$name)))
  # a$command <- "mv"
  # a$fastq <- paste0(a$`Comment[ENA_RUN]`,".fastq.gz")
  # a$newname <- gsub(".","_",a$name,fixed = TRUE)
  # a$newname <- paste0(a$newname,"_R1.fastq.gz")
  # lapply(a$newname,function(x)strsplit(x,split = ".",fix=TRUE)[[1]][2]) %>% unlist() %>% sort() %>% unique()
  
  count = fread(paste0("~/Pseudo/Data/Seqdata/Illumina/",Species,"/Kaessmann.Mergecount.raw.txt")) %>% as.data.frame()
  row.names(count) = count$V1
  count = count[,-1]
  
  name = colnames(count) %>% as.data.frame()
  name$change = a[match(name$.,a$`Comment[ENA_RUN]`),5]
  for (i in 1:nrow(a)) {
    colnames(count) %<>% gsub(a[i,1],a[i,5],.)
  }
  all(colnames(count)==name$change)
  fwrite(count,file = paste0("~/Pseudo/Data/Seqdata/Illumina/",Species,"/Kaessmann.Mergecount.changed.txt"),
         col.names = TRUE,row.names = TRUE, quote = FALSE,sep = "\t")
  rm(list = ls());gc();rm(list = ls())
}

#** Mouse ------------------------------------
{
  rm(list = ls());gc();rm(list = ls())
  Species = "Mouse"
  dir = paste0("~/Pseudo/Data/Seqdata/Illumina/",Species)
  a = grep("sdrf.txt$",list.files(dir),value=TRUE) %>% sapply(.,function(x){paste(dir,x,sep='/')}) %>%
    lapply(., function(x){ fread(x, header=T,sep = "\t")})
  a = a[[1]][,c(30,34)] %>% as.data.frame()
  a$V3 <- lapply(a$`Derived Array Data File`,function(x)strsplit(x,split = "Mouse.")[[1]][2]) %>% unlist()
  a$V4 <- lapply(a$V3,function(x)strsplit(x,split = ".sorted")[[1]][1]) %>% unlist()
  a <- a[order(a$V4),]
  f <- table(a$V4) %>% as.data.frame() 
  f$Var1 <- as.character(f$Var1)
  f <- f[order(f$Var1),]
  q <- c()
  for (i in 1:nrow(f)) {
    q <- c(q,paste0(f[i,1],"_rep",1:(f[i,2])))
  }
  a$name <- q
  a$name2 <- lapply(a$name,function(x)strsplit(x,split = "_rep")[[1]][1]) %>% unlist()
  table(a$V4==a$name2) %>% as.numeric() == length(unique(sort(a$name)))
  a$command <- "mv"
  a$fastq <- paste0(a$Comment.ENA_RUN.,".fastq.gz")
  a$newname <- a$name
  a$newname <- sub(".",".s",a$newname,fixed = TRUE)
  a$newname <- gsub(".5","",a$newname,fixed = TRUE)
  a$newname <- sub(".","_",a$newname,fixed = TRUE) #run two times
  a$newname <- sub(".","_",a$newname,fixed = TRUE)
  
  count = fread(paste0("~/Pseudo/Data/Seqdata/Illumina/",Species,"/Kaessmann.Mergecount.raw.txt")) %>% as.data.frame()
  row.names(count) = count$V1
  count = count[,-1]
  
  name = colnames(count) %>% as.data.frame()
  name$change = a[match(name$.,a$`Comment[ENA_RUN]`),9]
  for (i in 1:nrow(a)) {
    colnames(count) %<>% gsub(a[i,1],a[i,9],.)
  }
  all(colnames(count)==name$change)
  fwrite(count,file = paste0("~/Pseudo/Data/Seqdata/Illumina/",Species,"/Kaessmann.Mergecount.changed.txt"),
         col.names = TRUE,row.names = TRUE, quote = FALSE,sep = "\t")
  rm(list = ls());gc();rm(list = ls())
}


# Merge our own data ------------------------------------
rm(list = ls());gc();rm(list = ls())
dir = "/NAS/qians/NewCount/Mouse/Ourdata"
Files = grep("count$",list.files(dir),value=TRUE)
filePath <- sapply(Files,function(x){paste(dir,x,sep='/')})
data <- lapply(filePath, function(x){ fread(x, header=T,sep = "\t")})
a = data[[1]] %>% as.data.frame()
a = a[,c(1,ncol(a))]
colnames(a)[2] %<>% gsub(".q30.bam","",.)

for (i in 2:length(data)) {
  b <- data[[i]] %>% as.data.frame()
  b = b[,c(1,ncol(b))]
  colnames(b)[2] %<>% gsub(".q30.bam","",.)
  a <- merge(a,b,by="Geneid")
}
row.names(a) = a$Geneid
colnames(a) %<>% gsub("intestine","colon",.)
a = a[,-1]
fwrite(a,file = paste0("~/Pseudo/Data/Seqdata/Illumina/Mouse/Ourdata.Mergecount.txt"),
       col.names = TRUE,row.names = TRUE, quote = FALSE,sep = "\t")