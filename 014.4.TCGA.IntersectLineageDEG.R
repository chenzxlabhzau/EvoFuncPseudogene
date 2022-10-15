rm(list = ls());gc();rm(list = ls())
Num = "014.4."

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

# 1.All DEGs --------------------------------------------------------------
wd <- "~/Pseudo/Data/Seqdata/TCGA/DEG"   #current directory
directory = file.path(wd)
Files <- grep("DEG.csv$",list.files(directory),value=TRUE)
filePath <- sapply(Files, function(x){paste(directory,x,sep='/')})
data <- lapply(filePath, function(x){ fread(x)})
count = matrix(0,nrow = length(data), ncol = 3) %>% as.data.frame()
colnames(count) = c("Type","Number","Ratio") #pseudo number
#colnames(a) <- c("Gene",strsplit(colnames(a)[1],split=".",fix=TRUE)[[1]][1]) # geneID, sample name
for (i in 1:length(data)) {
  a <- data[i] %>% as.data.frame()
  count[i,] = c(colnames(a)[1],
                nrow(a[grepl("pseudogene",a[,8]),]),
                nrow(a[grepl("pseudogene",a[,8]),])/nrow(a)*100)
}
count$Type %<>% gsub("TCGA.","",.) %>% gsub(".DEG.csv.V1","",.) 
count$Number %<>% as.numeric()
count$Ratio %<>% as.numeric()
ggplot(count,aes(Type,Number))+geom_bar(stat = "identity",fill="#40a49a")+
  #geom_text(aes(label=Number), vjust=0.5 ,hjust=0.5) +
  xlab("Cancer type")+theme_classic()+
  theme(axis.text.x= element_text(size=12),axis.text.y= element_text(size=12),
        axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,2400),name = "Normalized iBAQ (log10)", 
                     sec.axis = sec_axis(trans = ~./24,
                                         name = "X:A ratio"))+coord_flip()+
    geom_point(aes(Type,Ratio*24))

p1 = ggplot(count,aes(Type,Number))+geom_bar(stat = "identity",fill="#40a49a")+
  geom_text(aes(label=Number), vjust=0.5 ,hjust=0.5) +
  theme_classic()+
  theme(axis.text.x= element_text(size=12),axis.text.y= element_blank(),
        axis.title.x = element_text(size=14), axis.title.y = element_blank(),axis.ticks.y = element_blank())+
  scale_y_continuous(expand = c(0, 0),limits = c(0,2400))+coord_flip()

p2 = ggplot(count,aes(Type,Ratio))+geom_point(size=3)+geom_line(aes(group=1))+
  xlab("Cancer types")+ylab("Percentage")+theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x= element_text(size=12),axis.text.y= element_text(size=12),
        axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,40))+coord_flip()
p3 = p2+p1+ plot_layout(ncol = 2, width = c(1, 1.5))
ggsave(p3,filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"PropNumber.DEGs.pdf")),
       device = "pdf",width = 5,height = 5)

# 2.upregulated DEGs --------------------------------------------------------------
count = matrix(0,nrow = length(data), ncol = 5) %>% as.data.frame()
colnames(count) = c("Type","Downpseu","Downall","Uppseu","Upall") #pseudo number
#colnames(a) <- c("Gene",strsplit(colnames(a)[1],split=".",fix=TRUE)[[1]][1]) # geneID, sample name
for (i in 1:length(data)) {
  a <- data[i] %>% as.data.frame()
  count[i,] = c(colnames(a)[1],
                #nrow(a[grepl("pseudogene",a[,8]),]),
                #nrow(a),
                nrow(a[grepl("pseudogene",a[,8]) & a[,3] <0,]),
                sum(a[,3]<0),
                nrow(a[grepl("pseudogene",a[,8]) & a[,3] >0,]),
                sum(a[,3]>0))
}
count$Type %<>% gsub("TCGA.","",.) %>% gsub(".DEG.csv.V1","",.) 
count[,-1] = apply(count[,-1],2, as.numeric)
table(count$Downpseu/count$Downall < count$Uppseu/count$Upall)
testP = c()
for (i in 1:nrow(count)) {
  matrix(as.numeric(count[i,2:5]),nrow = 2) %>% fisher.test()
  x = matrix(as.numeric(count[i,2:5]),nrow = 2) %>% fisher.test()
  testP[i] = x$p.value
}
count$Down = count$Downpseu / count$Downall * 100
count$Up = count$Uppseu / count$Upall * 100
row.names(count) = count$Type
count %<>% tidyr::pivot_longer(.,cols=6:7)
count$name %<>% factor(.,levels = rev(unique(.)))
ggplot(count,aes(Type,value))+geom_point(size=3)+geom_line(aes(group=name,linetype=name))+
  xlab("Cancer types")+ylab("Percentage")+theme_bw()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x= element_text(size=12,angle = 45,hjust = 1,vjust = 1),axis.text.y= element_text(size=12),
        axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))+
  geom_text(aes(x,y,label =lab),
            data = data.frame(x = 1:17,
                              y = ((testP>=0.05)/1)+60,
                              lab = cut(testP,breaks = c(1,0.05,0.01,0.001,0),
                                        labels = rev(c("n.s.","*","**","***")))))+
  theme(legend.title = element_blank(),legend.position="bottom")
ggsave(filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"Prop.DEGs.AllUp.pdf")),
         device = "pdf",width = 7,height = 5)
  
ggplot(count,aes(Type,value))+geom_point(size=3)+geom_line(aes(group=name,linetype=name))+
  xlab("Cancer types")+ylab("Percentage")+theme_bw()+
  theme(legend.key.size = unit(0.4, 'cm'),legend.position = c(0.65,0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x= element_text(size=12),axis.text.y= element_blank(),
        axis.title.x = element_text(size=14), axis.title.y = element_blank(),axis.ticks.y = element_blank())+
  geom_text(aes(x,y,label =lab),
            data = data.frame(x = 1:17,
                              y = ((testP>=0.05)/1)+56,
                              lab = cut(testP,breaks = c(1,0.05,0.01,0.001,0),
                                        labels = rev(c("n.s.","*","**","***")))))+
  coord_flip()+theme(legend.title = element_blank())
ggsave(filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"Prop.DEGs.AllUp.coordflip.pdf")),
       device = "pdf",width = 3,height = 5)

# 3.Intersect Lineage with DEG --------------------------------------------------------------
freq = matrix(0,nrow = length(data), ncol = 4) %>% as.data.frame()
colnames(freq) = c("Type","Number.of.Lineage","Number.of.DEG","Ratio")
lineage = fread("~/Pseudo/Data/Seqdata/TCGA/Expressbreadth/LineageSpcific.pseudo.csv") %>% as.data.frame()
lineage$dis %<>% lapply(.,function(x)strsplit(x,split=" ")[[1]][1]) %>% unlist()

for (i in 1:length(data)) {
  a <- data[i] %>% as.data.frame()
  freq[i,1:3] = c(colnames(a)[1],
                  nrow(a),
                  nrow(a[grepl("pseudogene",a[,8]),]))
}
freq$Type %<>% gsub("TCGA.","",.) %>% gsub(".DEG.csv.V1","",.) 


for (i in 1:length(data)) {
  a <- data[i] %>% as.data.frame()
  focal.cancer = colnames(a)[1] %>% gsub("TCGA.","",.) %>% gsub(".DEG.csv.V1","",.) 
  x = dplyr::filter(lineage, dis == focal.cancer)
  freq[i,1:4] = c(colnames(a)[1],length(x$V1),sum(x$V1 %in% a[,1]),sum(x$V1 %in% a[,1])/length(x$V1))
}
freq$Type %<>% gsub("TCGA.","",.) %>% gsub(".DEG.csv.V1","",.) 
freq$Ratio %<>% as.numeric()
write.table(freq,file = file.path("~/Pseudo/Result/TCGA/Savedata/",paste0(Num,"IntersectLineageDEG.csv")),
            quote = FALSE,sep = ",",row.names = FALSE)

# 4.Intersect Lineage with UP DEG --------------------------------------------------------------
freq = matrix(0,nrow = length(data), ncol = 4) %>% as.data.frame()
colnames(freq) = c("Type","Number.of.Lineage","Number.of.UP","Ratio")
lineage = fread("~/Pseudo/Data/Seqdata/TCGA/Expressbreadth/LineageSpcific.pseudo.csv") %>% as.data.frame()
lineage$dis %<>% lapply(.,function(x)strsplit(x,split=" ")[[1]][1]) %>% unlist()

for (i in 1:length(data)) {
  a <- data[i] %>% as.data.frame()
  focal.cancer = colnames(a)[1] %>% gsub("TCGA.","",.) %>% gsub(".DEG.csv.V1","",.) 
  x = dplyr::filter(lineage, dis == focal.cancer)
  freq[i,1:4] = c(colnames(a)[1],
                  length(x$V1),
                  sum(x$V1 %in% a[a[,3]>0,1]),
                  sum(x$V1 %in% a[a[,3]>0,1])/length(x$V1))
}
freq$Type %<>% gsub("TCGA.","",.) %>% gsub(".DEG.csv.V1","",.) 
freq$Ratio %<>% as.numeric()
write.table(freq,file = file.path("~/Pseudo/Result/TCGA/Savedata/",paste0(Num,"IntersectLineageDEG.csv")),
            quote = FALSE,sep = ",",row.names = FALSE)

