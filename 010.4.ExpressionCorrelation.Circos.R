rm(list = ls());gc();rm(list = ls())
Num = "010.4."
library(circlize)
#ref: https://jokergoo.github.io/circlize_book/book/genomic-introduction.html
Num = "010.4."
S = "Human"
pair = fread(file = file.path("~/Pseudo/Result",S,"Savedata","expressionPair.sig.csv"),nThread = 5) %>% 
  as.data.frame()
b = read.csv(paste0("~/Pseudo/Data/Ref/",S,"/","geneHuman.bed"),header = FALSE,sep = "\t")
gene = c(pair$pseu,pair$coding) 

circos.initializeWithIdeogram()
circos.initializeWithIdeogram(species = "hg38")
text(0, 0, "default", cex = 1)


set.seed(999)
bed = generateRandomBed(nr = 1000, nc = 1, fun = function(k) rnorm(k, 0, 0.5))
head(bed)
circos.par(track.height = 0.1,start.degree = 90,gap.degree = c(rep(2,20),2,2,2,10),canvas.xlim = c(-1,1.5))
circos.genomicInitialize(bed)
circos.genomicTrack(bed, panel.fun = function(region,value,...){
  circos.genomicPoints(region,value,numeric.column = 1,pch = 16,col = 'red',cex = 0.5)
})
circos.genomicTrack(bed, panel.fun = function(region,value,...){
  circos.genomicLines(region,value,numeric.column = 2,col = 'green')
})

circos.genomicTrack(bed, panel.fun = function(region,value,...){
  circos.genomicLines(region,value,col = 1:4)
})
col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))
circos.genomicTrack(bed, panel.fun = function(region,value,...){
  circos.genomicRect(region,value,col = col_fun(value[[1]]),border = NA)
})
circos.genomicTrack(bed, panel.fun = function(region,value,...){
  circos.genomicRect(region,value,col = "red",border = NA)
})


circos.track(ylim = c(0, 0.1), track.height = cm_h(0.1),
             track.margin = c(0, mm_h(0.1)),
             panel.fun = function(x, y) {
               xcenter = get.cell.meta.data("xcenter")
               circos.lines(c(xcenter, xcenter), c(0, cm_y(0.1)), col = "red")
             })
circos.clear()



col_meth = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
circos.heatmap(mat_meth, split = km, col = col_meth, track.height = 0.12)

col_direction = c("hyper" = "red", "hypo" = "blue")
circos.heatmap(direction, col = col_direction, track.height = 0.01)

col_expr = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
circos.heatmap(mat_expr, col = col_expr, track.height = 0.12)

col_pvalue = colorRamp2(c(0, 2, 4), c("white", "white", "red"))
circos.heatmap(cor_pvalue, col = col_pvalue, track.height = 0.01)

library(RColorBrewer)
col_gene_type = structure(brewer.pal(length(unique(gene_type)), "Set3"), names = unique(gene_type))
circos.heatmap(gene_type, col = col_gene_type, track.height = 0.01)

col_anno_gene = structure(brewer.pal(length(unique(anno_gene)), "Set1"), names = unique(anno_gene))
circos.heatmap(anno_gene, col = col_anno_gene, track.height = 0.01) 

col_dist = colorRamp2(c(0, 10000), c("black", "white"))
circos.heatmap(dist, col = col_dist, track.height = 0.01)

col_enhancer = colorRamp2(c(0, 1), c("white", "orange"))
circos.heatmap(anno_enhancer, col = col_enhancer, track.height = 0.3)




sectors = letters[1:10]
circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0))
circos.initialize(sectors, xlim = cbind(rep(0, 10), runif(10, 0.5, 1.5)))
circos.track(ylim = c(0, 1), track.height = mm_h(5),
             panel.fun = function(x, y) {
               circos.lines(c(0, 0 + mm_x(5)), c(0.5, 0.5), col = "blue")
             })
circos.track(ylim = c(0, 1), track.height = cm_h(1),
             track.margin = c(0, mm_h(2)),
             panel.fun = function(x, y) {
               xcenter = get.cell.meta.data("xcenter")
               circos.lines(c(xcenter, xcenter), c(0, cm_y(1)), col = "red")
             })
circos.clear()


#
S = "Human"
library(HelloRanges)
b = read.csv(paste0("~/Pseudo/Data/Ref/",S,"/","geneHuman.bed"),header = FALSE,sep = "\t")
b %<>% dplyr::filter(.,V1 %in% c(1:50,"X"))
b$V1 %<>% factor(.,levels = c(1:(length(unique(b$V1))-1),"X"))
circos.par(track.height = 0.1,start.degree = 82.5,
           gap.degree = c(rep(2,20),2,2,16),canvas.xlim = c(-1,1.5))
circos.initialize(factors = b$V1, x = b$V2)
circos.track(factors = b$V1,ylim = c(0,1),track.height =0.05,
             bg.col=c("#c77364","#ce8f5c","#7bac80","#75a5c0","#b5181a","#b72d9e",
                      "#e4cccb","#f6f0eb","#e8c755","#d5d456","#cfe74f","#39cf40",
                      "#3e6926","#0a0883","#6686c6","#a0c1db","#49ddd0","#e0f8f7",
                      "#651d91","#9d3eaf","#b9a3d0","#5b5a5b","#8f8f8f"),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(5), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })

pair = fread(file = file.path("~/Pseudo/Result",S,"Savedata","expressionPair.sig.csv"),nThread = 5) %>% 
  as.data.frame()
#for coding
cod = pair$coding %>% table() %>%  as.data.frame()
cod[,c(3:5)] = b[match(cod$.,b$V4),1:3]
cod = cod[,c(3:5,2)] %>% na.omit()
cod$V2 = cod$V2 - (cod$V3- cod$V2)*2
cod$V3 = cod$V3 + (cod$V3- cod$V2)
cod[cod$V2<0,"V2"]=0
fwrite(cod,file = file.path("~/Pseudo/Result",S,"Savedata","pair.coding.widen5.bed"),col.names = FALSE,sep = "\t")
#cd ~/Pseudo/Result/Human/Savedata
#awk -F "\t" '{print $1"\t""0""\t"$2}' ~/Pseudo/Data/Ref/Human/Human.length | bedtools subtract -a - -b ~/Pseudo/Result/Human/Savedata/pair.coding.widen5.bed > Human.inter4pairedcoding.bed
intercod = read.csv(file.path("~/Pseudo/Result",S,"Savedata",paste0(S,".inter4pairedcoding.bed")),header = FALSE,sep = "\t")
intercod$Freq = 0
cod = rbind(cod,intercod)
cod = cod[order(cod$V1,cod$V2),]
col_fun = colorRamp2(c(1, 5, 100), c("white", "black", "black"))
circos.genomicHeatmap(cod, numeric.column = 4, heatmap_height = 0.06, col = col_fun,connection_height = NULL)
#for pseudogene
pseu = pair$pseu %>% table() %>% as.data.frame()
pseu[,c(3:5)] = b[match(pseu$.,b$V4),1:3]
pseu = pseu[,c(3:5,2)] %>% na.omit()
pseu$V2 = pseu$V2 - (pseu$V3- pseu$V2)*2
pseu$V3 = pseu$V3 + (pseu$V3- pseu$V2)
pseu[pseu$V2<0,"V2"]=0
fwrite(pseu,file = file.path("~/Pseudo/Result",S,"Savedata","pair.pseudo.widen5.bed"),col.names = FALSE,sep = "\t")
#cd ~/Pseudo/Result/Human/Savedata
#awk -F "\t" '{print $1"\t""0""\t"$2}' ~/Pseudo/Data/Ref/Human/Human.length | bedtools subtract -a - -b ~/Pseudo/Result/Human/Savedata/pair.pseudo.widen5.bed > Human.inter4pairedpseudo.bed
interpseu = read.csv(file.path("~/Pseudo/Result",S,"Savedata",paste0(S,".inter4pairedpseudo.bed")),header = FALSE,sep = "\t")
interpseu$Freq = 0
pseu = rbind(pseu,interpseu)
pseu = pseu[order(pseu$V1,pseu$V2),]
col_fun = colorRamp2(c(1, 5, 100), c("white", "black", "black"))
circos.genomicHeatmap(pseu, numeric.column = 4, heatmap_height = 0.06, col = col_fun,connection_height = NULL)
#for link
pair[,paste0(c("c","s","e"),1)] = b[match(pair$coding,b$V4),1:3]
pair[,paste0(c("c","s","e"),2)] = b[match(pair$pseu,b$V4),1:3]
pair %<>% na.omit()
ddpseu = fread(paste0("~/Pseudo/Result/",S,"/Savedata/DDG/all.ddg.csv"),header = FALSE)
pair$type = "No"
pair[pair$pseu %in% ddpseu$V1,"type"] = "Yes"
bed1 = pair[,paste0(c("c","s","e"),1)]
bed2 = pair[,paste0(c("c","s","e"),2)]
circos.genomicLink(pair[pair$type=="No",paste0(c("c","s","e"),1)],
                   pair[pair$type=="No",paste0(c("c","s","e"),2)],
                   col = rand_color(1,transparency = 0.99))
circos.genomicLink(pair[pair$type=="Yes",paste0(c("c","s","e"),1)],
                   pair[pair$type=="Yes",paste0(c("c","s","e"),2)],
                   col = rand_color(1,transparency = 0.99))
circos.clear()
############################# final
pdf(file =  paste0("/home/qians/Pseudo/Result/",S,"/Picture/",Num,S,".coexpressedpair.circos.pdf"),width = 6,height = 5.5)
circos.par(track.height = 0.1,start.degree = 82.5,
           gap.degree = c(rep(2,20),2,2,16),canvas.xlim = c(-1,1.5))
circos.initialize(factors = b$V1, x = b$V2)
circos.track(factors = b$V1,ylim = c(0,1),track.height =0.05,
             bg.col=c("#c77364","#ce8f5c","#7bac80","#75a5c0","#b5181a","#b72d9e",
                      "#e4cccb","#f6f0eb","#e8c755","#d5d456","#cfe74f","#39cf40",
                      "#3e6926","#0a0883","#6686c6","#a0c1db","#49ddd0","#e0f8f7",
                      "#651d91","#9d3eaf","#b9a3d0","#5b5a5b","#8f8f8f"),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(5), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
circos.genomicHeatmap(cod, numeric.column = 4, heatmap_height = 0.06, col = col_fun,connection_height = NULL)
circos.genomicHeatmap(pseu, numeric.column = 4, heatmap_height = 0.06, col = col_fun,connection_height = NULL)
circos.genomicLink(pair[pair$type=="No",paste0(c("c","s","e"),1)],
                   pair[pair$type=="No",paste0(c("c","s","e"),2)],
                   col = rand_color(1,transparency = 0.98))
circos.genomicLink(pair[pair$type=="Yes",paste0(c("c","s","e"),1)],
                   pair[pair$type=="Yes",paste0(c("c","s","e"),2)],
                   col = rand_color(1,transparency = 0.98))
circos.clear()
dev.off()
#############################

bed1 = generateRandomBed(nr = 100)
bed1 = bed1[sample(nrow(bed1), 20), ]
bed2 = generateRandomBed(nr = 100)
bed2 = bed2[sample(nrow(bed2), 20), ]

circos.initializeWithIdeogram()
circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.99), 
                   border = NA)




pseu = pseu[,c(3:5,2)] %>% na.omit()
col_fun = colorRamp2(c(2, 5, 100), c("white", "black", "red"))
circos.genomicHeatmap(pseu, numeric.column = 4, heatmap_height = 0.06, col = col_fun,connection_height = NULL)


circos.genomicTrack(cod, numeric.column = 4, panel.fun = function(region,value,...){
  circos.genomicRect(region,value,col = col_fun(value[[1]]),border = NA)
})



circos.genomicTrack(cod, panel.fun = function(region,value,...){
  circos.genomicHeatmap(cod, numeric.column = 4, heatmap_height = 0.06, col = col_fun,connection_height = NULL)
})

circos.genomicTrack(cod, panel.fun = function(region,value,...){
  circos.dendrogram(region,value)
})


col_fun = colorRamp2(c(-1, 0, 1), c("#d7191c", "#f7f7f7", "#2c7bb6"))
circos.genomicTrack(cod, numeric.column = 4, 
                    panel.fun = function(region,value,...){
                      circos.genomicRect(region, value, col = col_fun(value[[1]]),
                                         #ybottom = 0,ytop = 0.1,
                                         border = NA,ylim = CELL_META$ylim)
                    })

circos.heatmap(anno_enhancer, col = col_enhancer, track.height = 0.3)

set.seed(123)
mat1 = rbind(cbind(matrix(rnorm(50*5, mean = 1), nr = 50), 
                   matrix(rnorm(50*5, mean = -1), nr = 50)),
             cbind(matrix(rnorm(50*5, mean = -1), nr = 50), 
                   matrix(rnorm(50*5, mean = 1), nr = 50))
)
rownames(mat1) = paste0("R", 1:100)
colnames(mat1) = paste0("C", 1:10)
mat1 = mat1[sample(100, 100), ] # randomly permute rows
split = sample(letters[1:5], 100, replace = TRUE)
spilt = factor(split, levels = letters[1:5])
col_fun1 = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
circos.heatmap(mat1, split = split, col = col_fun1)
circos.clear()



circos.clear()
#
set.seed(123)
circos.initializeWithIdeogram(plotType = NULL,start.degree = 90,)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)



if (TRUE) {
  set.seed(999)
  n = 1000
  df = data.frame(factors = sample(letters[1:3], n, replace = TRUE),
                  x = abs(rnorm(n)*100), y = runif(n))
  circos.par("track.height" = 0.1)
  circos.par(track.height = 0.1,start.degree = 90,gap.degree = c(2,2,10),canvas.xlim = c(-1,1.5))
  circos.initialize(factors = df$factors, x = df$x)
  
  circos.track(factors = df$factors, y = df$y,bg.col=c("#957244","#FFDD89","#F26223"),
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(5), 
                             CELL_META$sector.index)
                 circos.axis(labels.cex = 0.6)
               })
  col = rep(c("#FF0000", "#00FF00"), 4)
  circos.trackPoints(df$factors, df$x, df$y, col = col, pch = 16, cex = 0.5)
  circos.text(-1, 0.5, "text", sector.index = "a", track.index = 1)
  
  bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
  circos.trackHist(df$factors, df$x, bin.size = 0.2, bg.col = bgcol, col = NA)
  
  circos.track(factors = df$factors, x = df$x, y = df$y,
               panel.fun = function(x, y) {
                 ind = sample(length(x), 10)
                 x2 = x[ind]
                 y2 = y[ind]
                 od = order(x2)
                 circos.lines(x2[od], y2[od])
               })
  
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 0.1)
    n_breaks = length(breaks)
    circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
                breaks[-1], rep(ylim[2], n_breaks - 1),
                col = rand_color(n_breaks), border = NA)
  })
  
  circos.link("a", 0, "b", 0, h = 0.4)
  circos.link("c", c(-0.5, 0.5), "d", c(-0.5,0.5), col = "red",
              border = "blue", h = 0.2)
  circos.link("e", 0, "g", c(-1,1), col = "green", border = "black", lwd = 2, lty = 2)
}




