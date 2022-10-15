rm(list = ls());gc();rm(list = ls())
Num = "014.5."

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

wd <- "~/Pseudo/Data/Seqdata/TCGA/DEG"   #current directory
directory = file.path(wd)
Files <- grep("DEG.csv$",list.files(directory),value=TRUE)
filePath <- sapply(Files, function(x){paste(directory,x,sep='/')})
data <- lapply(filePath, function(x){ fread(x)})
count = matrix(0,nrow = length(data), ncol = 9) %>% as.data.frame()
colnames(count) = c("Type","DEpseu","Allpseu","DynamicDE","Alldynamic","CI1","CI2","OR","Pval") 
#colnames(a) <- c("Gene",strsplit(colnames(a)[1],split=".",fix=TRUE)[[1]][1]) # geneID, sample name
for (i in 1:length(data)) {
  a <- data[i] %>% as.data.frame()
  count[i,1:2] = c(colnames(a)[1],
                nrow(a[grepl("pseudogene",a[,8]),]))
}
count$Type %<>% gsub("TCGA.","",.) %>% gsub(".DEG.csv.V1","",.) 
count$DEpseu %<>% as.numeric()
count$Allpseu = 15244

ddg = fread("~/Pseudo/Result/Human/Savedata/DDG/all.ddg.csv",header = FALSE) %>% as.data.frame()
for (i in 1:length(data)) {
  a <- data[i] %>% as.data.frame()
  count[i,4] = sum(a[,1] %in% ddg$V1)
}

count$Alldynamic = length(unique(ddg$V1))

count$All = count$DEpseu / count$Allpseu * 100
count$Dynamic = count$DynamicDE / count$Alldynamic * 100
row.names(count) = count$Type
count = count[rev(1:nrow(count)),]
library(ggplotify)
p1 = ComplexHeatmap::pheatmap(count[,10:11],#annotation_legend = T,
                         border_color = 'black',#设定每个格子边框的颜色，border=F则无边框
                         cluster_rows = F, #对行聚类
                         cluster_cols = F,
                         #column_names_rot = 45,
                         #column_title = "Proportion of pseudogenes",
                         row_title = "Cancer type") %>% as.ggplot
ggsave(filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"ProporDEG.AllDynamic.pdf")),
       p1,device = "pdf",width = 2.5,height = 4)

fisher.test(matrix(as.numeric(count[1,2:5]),nrow = 2))

for (i in 1:length(data)) {
  x = fisher.test(matrix(as.numeric(count[i,c(4,5,2,3)]),nrow = 2))
  count[i,"CI1"] = x$conf.int[1]
  count[i,"CI2"] = x$conf.int[2]
  count[i,"OR"] = x$estimate
  count[i,"Pval"] = x$p.value
}

p2 = ggplot(count,aes(Type,OR))+geom_point()+
  geom_errorbar(aes(x = Type, ymax=CI1, ymin=CI2,width =0.3),position = position_dodge(0.5))+
  coord_flip()+theme_bw()+
  ylab("Odds ratio")+
  theme(panel.grid.major =element_blank(),#text=element_text(family="Times New Roman"),
        #title = element_text(family="Times New Roman"),
        #panel.grid.minor = element_blank(),
        axis.text.x= element_text(size=12, color="black"),axis.text.y= element_blank(),
        axis.title.x = element_text(size=14), axis.title.y = element_blank(),axis.ticks.y = element_blank())+
  geom_hline(yintercept = 1,color="#7aaed8")
ggsave(filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"OddsRatio.AllDynamic.pdf")),
       p2,device = "pdf",width = 2.5,height = 4)

pdf(file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"test.pdf")),width = 5,height = 5.5)
p1
p2
dev.off()
