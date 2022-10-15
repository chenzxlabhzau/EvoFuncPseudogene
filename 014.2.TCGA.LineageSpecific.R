rm(list = ls());gc();rm(list = ls())
Num = "014.2."

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

i = 1
focal.cancer = paste0("TCGA-",all[i])
a = fread(paste0("~/Pseudo/Data/Seqdata/TCGA/MergeCount/",focal.cancer,".mergecount.csv")) %>% as.data.frame()
row.names(a) = a$V1
a = a[,-1]
a %<>% dplyr::select(.,ends_with("TP"))
cpm = apply(a, 2, function(x)x/sum(x) * 10^6) %>% as.data.frame()
cpm$mean = apply(cpm,1,mean)
mtx = cpm[,"mean",drop=FALSE]
colnames(mtx) = paste0(focal.cancer,"(",ncol(a),")")
rm(a);rm(cpm)

for (i in 2:length(all)) {
  focal.cancer = paste0("TCGA-",all[i])
  a = fread(paste0("~/Pseudo/Data/Seqdata/TCGA/MergeCount/",focal.cancer,".mergecount.csv")) %>% as.data.frame()
  row.names(a) = a$V1
  a = a[,-1]
  a %<>% dplyr::select(.,ends_with("TP"))
  cpm = apply(a, 2, function(x)x/sum(x) * 10^6) %>% as.data.frame()
  cpm$mean = apply(cpm,1,mean)
  cpm = cpm[,"mean",drop=FALSE]
  #colnames(cpm) = paste0(focal.cancer,"(",ncol(a),")")
  mtx[,paste0(focal.cancer,"(",ncol(a),")")] = cpm[match(row.names(mtx),row.names(cpm)),1]
}
colnames(mtx) %<>% gsub("TCGA-","",.) %>% gsub("("," (",.,fixed = TRUE)

#1.use tau for filter lineage-specific genes
if (FALSE) {
  #tau <- function(x){
  #  x <- x[apply(x, 1, function(y)sum(y >0.1 ) >=1),]
  #  x <-  as.data.frame(t(apply(x, 1,function(y){1-y/max(y)})))
  #  x$tau <- apply(x,1,function(y)sum(y)/(length(y)-1))
  #  return(as.data.frame(x))
  #}
  tmp <- tau(mtx)
  x = mtx
  x = x[apply(x, 1, sum)>0,]
  x$tau = tmp[match(row.names(x),row.names(tmp)),"tau"]
  x = x[x$tau>0.8,] %>% na.omit()
  pheatmap(x[,-ncol(x)],scale = "row",
           colorRampPalette(c("#f4f4c8", "#fd0c06"))(50),
           #annotation_colors = ann_colors,
           #annotation_row = row,
           annotation_legend = T,
           border_color = 'black',#设定每个格子边框的颜色，border=F则无边框
           cluster_rows = F, #对行聚类
           cluster_cols = F, #队列聚类
           show_colnames = T, #是否显示列名
           show_rownames = F, #是否显示行名
           main = "TCGA"
  )
}

#2.use expression level
lineage = mtx
lineage = lineage[apply(lineage, 1, sum)>0,]
lineage = apply(lineage, 1, function(x)x/sum(x)*100) %>% t() %>% as.data.frame()
lineage = lineage[apply(lineage, 1, max) >15,]
lineage = lineage[apply(lineage, 1, function(x)sort(x,decreasing = TRUE)[2]) < 5,]
lineage$dis = colnames(lineage)[apply(lineage,1, which.max)]
lineage$max = apply(lineage[,-ncol(lineage)],1,max)
lineage = lineage[order(lineage$dis,-lineage$max),]
b = read.csv("~/Pseudo/Data/Ref/Human/geneHuman.bed",header = FALSE,sep = "\t")
write.table(lineage[intersect(row.names(lineage) ,b[grepl("pseu",b$V5),4]),],
            file = "~/Pseudo/Data/Seqdata/TCGA/Expressbreadth/LineageSpcific.pseudo.csv",
            quote = FALSE,sep = "\t")
write.table(lineage,
            file = "~/Pseudo/Data/Seqdata/TCGA/Expressbreadth/LineageSpcific.all.csv",
            quote = FALSE,sep = "\t")
#colnames(lineage) %<>% gsub("TCGA-","",.) %>% gsub("("," (",.,fixed = TRUE)

# (All, Pseudo) + (Lineage, Ubiquitous) -----------------------------------
#* All + Lineage -----------------------------------
p1 = pheatmap(lineage[,-c(ncol(lineage),ncol(lineage)-1)],#scale = "row",
         colorRampPalette(c("#f4f4c8", "#ae1e1e"))(50),
         #annotation_colors = ann_colors,
         #annotation_row = row,
         annotation_legend = T,
         border_color = 'black',#设定每个格子边框的颜色，border=F则无边框
         cluster_rows = F, #对行聚类
         cluster_cols = F, #队列聚类
         show_colnames = T, #是否显示列名
         show_rownames = F, #是否显示行名
         main = "All genes"
)
ggsave(p1,filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"Allgene.LineageSpecific.pdf")),
       device = "pdf",width = 5,height = 4)
#* Pseudo + Lineage -----------------------------------
p2 = pheatmap(lineage[intersect(row.names(lineage) ,b[grepl("pseu",b$V5),4]),-c(ncol(lineage),ncol(lineage)-1)],#scale = "row",
              colorRampPalette(c("#f4f4c8", "#ae1e1e"))(50),
              #annotation_colors = ann_colors,
              #annotation_row = row,
              annotation_legend = T,
              border_color = 'black',#设定每个格子边框的颜色，border=F则无边框
              cluster_rows = F, #对行聚类
              cluster_cols = F, #对列聚类
              show_colnames = T, #是否显示列名
              show_rownames = F #是否显示行名
)
ggsave(p2,filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"Pseudogene.LineageSpecific.pdf")),
       device = "pdf",width = 5,height = 4)

#* All + Ubiquitous -----------------------------------
ubiquitous = mtx
ubiquitous = ubiquitous[apply(ubiquitous, 1, sum)>0,]
ubiquitous = apply(ubiquitous, 1, function(x)x/sum(x)*100) %>% t() %>% as.data.frame()
ubiquitous = ubiquitous[apply(ubiquitous, 1, max) <= 30,]
ubiquitous = ubiquitous[apply(ubiquitous, 1, function(x)sort(x,decreasing = TRUE)[5]) > 5,]
write.table(ubiquitous[intersect(row.names(ubiquitous) ,b[grepl("pseu",b$V5),4]),],
            file = "~/Pseudo/Data/Seqdata/TCGA/Expressbreadth/Ubiquitous.pseudo.csv",
            quote = FALSE,sep = "\t")
write.table(lineage,
            file = "~/Pseudo/Data/Seqdata/TCGA/Expressbreadth/Ubiquitous.all.csv",
            quote = FALSE,sep = "\t")
#ubiquitous = ubiquitous[apply(ubiquitous, 1, function(x)sort(x,decreasing = TRUE)[10]) > 3,]
#ubiquitous = ubiquitous[apply(ubiquitous, 1, function(x)sort(x,decreasing = TRUE)[2]) < 5,]
#ubiquitous$dis = colnames(ubiquitous)[apply(ubiquitous,1, which.max)]
#ubiquitous$max = apply(ubiquitous[,-ncol(ubiquitous)],1,max)
#ubiquitous = ubiquitous[order(ubiquitous$dis,-ubiquitous$max),]
p3 = pheatmap(ubiquitous,#scale = "row",
         colorRampPalette(c("#f4f4c8", "#fd0c06"))(50),
         #annotation_colors = ann_colors,
         #annotation_row = row,
         annotation_legend = T,
         border_color = 'black',#设定每个格子边框的颜色，border=F则无边框
         cluster_rows = F, #对行聚类
         cluster_cols = F, #队列聚类
         show_colnames = T, #是否显示列名
         show_rownames = F, #是否显示行名
         main = "All gene"
)
ggsave(p3,filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"All.Ubiquitous.pdf")),
       device = "pdf",width = 5,height = 4)
#* Pseudogene + Ubiquitous -----------------------------------
p4 = pheatmap(ubiquitous[intersect(row.names(ubiquitous) ,b[grepl("pseu",b$V5),4]),-c(ncol(ubiquitous),ncol(ubiquitous)-1)],#scale = "row",
              colorRampPalette(c("#f4f4c8", "#ae1e1e"))(50),
              #annotation_colors = ann_colors,
              #annotation_row = row,
              annotation_legend = T,
              border_color = 'black',#设定每个格子边框的颜色，border=F则无边框
              cluster_rows = F, #对行聚类
              cluster_cols = F, #队列聚类
              show_colnames = T, #是否显示列名
              show_rownames = F #是否显示行名
)
ggsave(p4,filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"Pseudogene.Ubiquitous.pdf")),
       device = "pdf",width = 5,height = 4)

# Intersect with DDG -----------------------------------
ddg = fread("~/Pseudo/Result/Human/Savedata/DDG/all.ddg.csv",header = FALSE) %>% as.data.frame()
#for lineage-specific pseudo
fisher.test(matrix(c(sum(row.names(lineage) %in% ddg$V1), #N of ddg
                   sum(row.names(lineage) %in% b[grepl("pseu",b$V5),4]), #N of pseudo
                   length(unique(ddg$V1)), #all ddg
                   length(b[grepl("pseu",b$V5),4])), #all pseudo
                   nrow = 2)) #p = 9e-04
#for ubiquitous pseudo
fisher.test(matrix(c(sum(row.names(ubiquitous) %in% ddg$V1), #N of ddg
                     sum(row.names(ubiquitous) %in% b[grepl("pseu",b$V5),4]), #N of pseudo
                     length(unique(ddg$V1)), #all ddg
                     length(b[grepl("pseu",b$V5),4])), #all pseudo
                   nrow = 2))
fish = data.frame("Specific"=sum(row.names(lineage) %in% ddg$V1)/
           sum(row.names(lineage) %in% b[grepl("pseu",b$V5),4])*100,
           "Ubiquitous"=sum(row.names(ubiquitous) %in% ddg$V1)/
             sum(row.names(ubiquitous) %in% b[grepl("pseu",b$V5),4])*100) %>% t() %>% as.data.frame()
fish$type = row.names(fish)

p5 = ggplot(fish,aes(type,V1,fill=type))+geom_bar(stat = "identity")+
  scale_fill_manual(values=c("#ae1e1e","#D08671"))+
  theme_classic()+ylab("(%) DDPs")+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=12),axis.title.x = element_blank(),
        axis.text.y = element_text(size=12),axis.title.y = element_text(size=14))+
  guides(fill=FALSE)+
  scale_y_continuous(expand = c(0, 0),limits = c(0,40))+
  #geom_hline(yintercept = length(unique(ddg$V1))/length(b[grepl("pseu",b$V5),4])*100,color="black",linetype = "dashed",size=0.8)+
  geom_segment(aes(x=1, y=37.5, xend=2, yend=37.5))+
  annotate("text", x=1.5, y=38, label="**",size=7) #p=0.004
ggsave(p5,filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"PropDynamic.PseudoType.pdf")),
       device = "pdf",width = 3,height = 4.5)

fisher.test(matrix(c(sum(row.names(lineage) %in% ddg$V1), #N of ddg
                     sum(row.names(lineage) %in% b[grepl("pseu",b$V5),4]), #N of pseudo
                     sum(row.names(ubiquitous) %in% ddg$V1), #N of ddg
                     sum(row.names(ubiquitous) %in% b[grepl("pseu",b$V5),4])), #all pseudo
                     nrow = 2)) #p=0.004

