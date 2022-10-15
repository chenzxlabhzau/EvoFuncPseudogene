rm(list = ls());gc();rm(list = ls())
Num = "004.5."

# Gene traffic ----
species = "Human"
a = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),header = TRUE,sep = ",",stringsAsFactors = FALSE)
pair = read.csv("~/Pseudo/Data/Ref/Human/pseudo.coding.pair.txt",header = TRUE,sep = "\t",stringsAsFactors = FALSE)
pair[,c("pseudochr","age")] = a[match(pair$Pseudogene,a$Human),c("chr","age")]
pair %<>% na.omit()

b = read.csv("~/Pseudo/Data/Ref/Human/geneHuman.bed",header = FALSE,sep = "\t",stringsAsFactors = FALSE)
b %<>% dplyr::filter(.,V1 %in% c(1:22,"X"))
pair$codingchr = b[match(pair$Parentgene,b$V4),1]
pair %<>% dplyr::filter(.,!codingchr %in% c("MT","Y")) %>% na.omit()
#pair[pair$codingchr!="X","codingchr"] = "A"
#pair[pair$pseudochr!="X","pseudochr"] = "A"
table(pair$pseudochr=="X",pair$codingchr=="X") %>% fisher.test()

# *For all psuedogene ----
# ** Parent coding ~ chr coding gene number ----
coding = table(pair$codingchr) %>% as.data.frame()
chrnumber = table(b[b$V5==" protein_coding",1]) %>% as.data.frame()
coding$chrnumber = chrnumber[match(coding$Var1,chrnumber$Var1),2]
coding$group = c(rep("A",22),"X")
ggplot(coding,aes(chrnumber,Freq,color=group))+geom_point()+geom_smooth(method = "lm") #Negative

# ** Pseudo ~ chr length ----
pseudo = table(pair$pseudochr) %>% as.data.frame()
chrlength = read.csv("~/Pseudo/Data/Ref/Human/Human.length",header = FALSE,sep = "\t")
pseudo$chrlength = chrlength[match(pseudo$Var1,chrlength$V1),2]
pseudo$group = c(rep("A",22),"X")
ggplot(pseudo,aes(chrlength,Freq,color=group))+geom_point()+geom_smooth(method = "lm") #No sig
 
# ** Expect/Observe ----
##ref:Extensive Gene Traffic on the Mammalian X Chromosome
#All pseu
##ratio
rlength = chrlength[chrlength$V1=="X",2]/ sum(chrlength[chrlength$V1 %in% c(1:22,"X"),2])
rgene = chrnumber[chrnumber$Var1=="X",2]/ sum(chrnumber[chrnumber$Var1 %in% c(1:22,"X"),2])
matrix(data = c(nrow(pair)*rgene,
                nrow(pair)*(1-rgene),
                sum(pair$codingchr=="X"),
                sum(pair$codingchr!="X")),
       nrow = 2) %>% fisher.test()

matrix(data = c(nrow(pair)*rlength,
                nrow(pair)*(1-rlength),
                sum(pair$pseudochr=="X"),
                sum(pair$pseudochr!="X")),
       nrow = 2) %>% fisher.test()
#Processed pseudo
pair$type = b[match(pair$Pseudogene,b$V4),5]
dim(pair[pair$type %in% c(" processed_pseudogene"," transcribed_processed_pseudogene"),])
matrix(data = c(nrow(pair[pair$type %in% c(" processed_pseudogene"," transcribed_processed_pseudogene"),])*rgene,
                nrow(pair[pair$type %in% c(" processed_pseudogene"," transcribed_processed_pseudogene"),])*(1-rgene),
                sum(pair[pair$type %in% c(" processed_pseudogene"," transcribed_processed_pseudogene"),"codingchr"]=="X"),
                sum(pair[pair$type %in% c(" processed_pseudogene"," transcribed_processed_pseudogene"),"codingchr"]!="X")),
       nrow = 2) %>% fisher.test()

matrix(data = c(nrow(pair[pair$type %in% c(" processed_pseudogene"," transcribed_processed_pseudogene"),])*rlength,
                nrow(pair[pair$type %in% c(" processed_pseudogene"," transcribed_processed_pseudogene"),])*(1-rlength),
                sum(pair[pair$type %in% c(" processed_pseudogene"," transcribed_processed_pseudogene"),"pseudochr"]=="X"),
                sum(pair[pair$type %in% c(" processed_pseudogene"," transcribed_processed_pseudogene"),"pseudochr"]!="X")),
       nrow = 2) %>% fisher.test()

# ** Traffic visual ----
## ALL
pair$pseudochr %<>% factor(.,levels = c(1:22,"X"))
pair$codingchr %<>% factor(.,levels = c(1:22,"X"))
m <- matrix(as.numeric(table(pair$pseudochr,pair$codingchr)),
            byrow = FALSE,
            nrow = 23, ncol = 23)
haircolors <- sort(unique(pair$pseudochr))
dimnames(m) <- list(have = haircolors,
                    prefer = haircolors)

groupColors <- c(rep("#8db8d5",22),"#9c0b1f")
#pdf(file =  paste0("/home/qians/Pseudo/Result/Human/Picture/",Num,"Genetraffic.pdf"),width = 6,height = 6)
chorddiag(m, groupColors = groupColors, groupnamePadding = 25,ticklabelFontsize = 8)

# ** Express bias ----
wd = "~/MamDC/Data/Seqdata/CardosoKaessmann2019NatureRNAseq/"
exp = fread(paste0(wd,species,"/",species,".allgene.fpkm.txt")) %>% as.data.frame()
colnames(exp) %<>% gsub("KidneyTestis_4w_Male","Kidney_4w_Male",.) %>%
  gsub("Forebrain","Brain",.) %>% gsub("Hindbrain","Cerebellum",.)
b$V5 %<>% gsub(" ","",.)
exp[,"type"] <- b[match(exp$V1,b$V4),5]
#exp %<>% dplyr::filter(.,grepl("pseudogene",type,ignore.case = T))
row.names(exp) = exp$V1
exp = exp[,-1]
exp = exp[apply(exp[,-ncol(exp)], 1, max)>0,] #filter 0FPKM in all tissues
exp[1:3,1:3]
exp$sample = colnames(exp)[apply(exp[,-ncol(exp)], 1, which.max)]
exp$tissue = lapply(exp$sample,function(x)(strsplit(x,"_")[[1]][1])) %>% unlist()

pair$tparent = exp[match(pair$Parentgene,row.names(exp)),"tissue"]
pair$tpseu = exp[match(pair$Pseudogene,row.names(exp)),"tissue"]
pair %<>% na.omit()
df = table(pair$tparent,pair$tpseu,paste0(pair$pseudochr=="X",pair$codingchr=="X"))
df[,,1][7,7]/sum(df[,,1][7,])
sum(df[,,1][1:6,7])/sum(df[,,1][7,])

  table(paste0(pair$codingchr=="X",pair$pseudochr=="X"),paste(pair$tparent,pair$tpseu,sep="2")) %>%
  as.data.frame()
df %<>% pivot_wider(.,names_from = "dis1",values_from="Freq")
df$dis1 = lapply(df$Var2,function(x)strsplit(as.character(x),split = "2",fixed = TRUE)[[1]][1]) %>% unlist()
