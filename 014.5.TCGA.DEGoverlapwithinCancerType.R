rm(list = ls());gc();rm(list = ls())
Num = "014.5."

wd <- "~/Pseudo/Data/Seqdata/TCGA/DEG"   #current directory
directory = file.path(wd)
Files <- grep("DEG.csv$",list.files(directory),value=TRUE)
filePath <- sapply(Files, function(x){paste(directory,x,sep='/')})
data <- lapply(filePath, function(x){ fread(x)})
names  = c()
for (i in 1:length(data)) {
  a <- data[i] %>% as.data.frame()
  a = a[grepl("pseu",a[,8]),]
  names = c(names, a[,1])
}
names %<>% unique()

mtx = matrix(0,nrow = length(names), ncol = length(data)) %>% as.data.frame()
row.names(mtx) = names
for (i in 1:length(data)) {
  a <- data[i] %>% as.data.frame()
  colnames(mtx)[i] = colnames(a)[1] %>% gsub("TCGA.","",.) %>% gsub(".DEG.csv.V1","",.) 
  #mtx[row.names(mtx) %in% a[,1], i] =1
  mtx[,i] = a[match(row.names(mtx),a[,1]),3]
}
mtx[is.na(mtx)] = 0

m = cor(mtx)
corrplot::corrplot(m)

# Dissimilarity matrix
d <- dist(as.data.frame(t(mtx)))

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

pheatmap(mtx, scale = "row",
         cluster_cols = T, #对列聚类
         cluster_rows = F,
         show_colnames = T, #是否显示列名
         show_rownames = F )#是否显示行名)

#Ref: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
library(corrplot)
pdf(file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"Corrplot.DEGs.pdf")),width = 5,height = 5.5)
corrplot(m, method = "color", col = rev(col(200)),
         type = "lower", order = "hclust", 
         number.cex = .7,
         #addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black",tl.srt = 90, # Text label color and rotation
         # Combine with significance
         sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)
dev.off()
