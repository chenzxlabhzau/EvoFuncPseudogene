rm(list = ls());gc();rm(list = ls())
Num = "015.6."

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
count = matrix(0,nrow = length(all),ncol = 30) %>% as.data.frame()
row.names(count) = all
colnames(count) =  paste(rep(c("B","CD4","CD8","Neutrophil","Macrophage","Dendritic"),each=5),
                      c("Freq","OR1","OR2","OR3","P"),sep = ".")

for (i in 1:length(all)) {
  focal.cancer = paste0("TCGA-",all[i]) 
  immu = fread(file =  paste0("/home/qians/Pseudo/Data/Seqdata/TCGA/ImmunoPseudo/",
                           focal.cancer,".ImmuPseudo.csv")) %>% 
    as.data.frame() %>% na.omit()
  #filter
  immu %<>% dplyr::filter(.,abs(sigValue) > 0.995 & padj < 0.05)
  
  infil = fread(file = paste0("~/Pseudo/Data/Seqdata/TCGA/Infiltration/",focal.cancer,".infiltration.cor.csv")) %>%
    as.data.frame() %>% na.omit()
  row.names(infil) = infil$V1
  infil = infil[,-1]
  
  for (j in 1:6) {
    #Freq
    count[i,5*j-4] = sum(row.names(infil)[abs(infil[,2*j-1])>0.3 & infil[,2*j]<0.05] %in% immu$lncRNA) / 
      sum(abs(infil[,2*j-1])>0.3 & infil[,2*j]<0.05)
    #Odds ratio
    x = matrix(c(sum(row.names(infil)[abs(infil[,2*j-1])>0.3 & infil[,2*j]<0.05] %in% immu$lncRNA),
             sum(abs(infil[,2*j-1])>0.3 & infil[,2*j]<0.05),
             length(unique(immu$lncRNA)),
             nrow(infil)),nrow = 2) %>% fisher.test()
    count[i,c(5*j-3,5*j-2,5*j-1,5*j)] = c(x$conf.int[1],x$estimate,x$conf.int[2],x$p.value)
  }
}


df = matrix(0,nrow = nrow(count)*6, ncol = 5+1) %>% as.data.frame()
colnames(df) = c("Freq","OR1","OR2","OR3","P","Type")
df[1:nrow(count),1:5] = as.matrix(count[,1:5])
df[1:nrow(count),6] = "B"
df[(nrow(count)+1):(nrow(count)*2),1:5] = as.matrix(count[,6:10])
df[(nrow(count)+1):(nrow(count)*2),6] = "CD4"
df[(nrow(count)*2+1):(nrow(count)*3),1:5] = as.matrix(count[,11:15])
df[(nrow(count)*2+1):(nrow(count)*3),6] = "CD8"
df[(nrow(count)*3+1):(nrow(count)*4),1:5] = as.matrix(count[,16:20])
df[(nrow(count)*3+1):(nrow(count)*4),6] = "Neutrophil"
df[(nrow(count)*4+1):(nrow(count)*5),1:5] = as.matrix(count[,21:25])
df[(nrow(count)*4+1):(nrow(count)*5),6] = "Macrophage"
df[(nrow(count)*5+1):(nrow(count)*6),1:5] = as.matrix(count[,26:30])
df[(nrow(count)*5+1):(nrow(count)*6),6] = "Dendritic"
df$Cancer = rep(row.names(count),6)
  
df[df$Cancer=="OV"&df$Type=="CD8",2:4]=0
df[df$Cancer=="BLCA"&df$Type=="B",2:4]=0
df$Type %<>% gsub("B","B cell",.) %>% gsub("CD4","CD4 cell",.)%>% gsub("CD8","CD8 cell",.)
ggplot(df,aes(Cancer,OR2))+geom_point()+facet_wrap(.~Type,scales = "free_x",ncol = 3)+
  geom_errorbar(aes(x = Cancer, ymax=OR1, ymin=OR3,width =0.01),position = position_dodge(0.5))+
  coord_flip()+theme_bw()+
  ylab("Odds ratio")+
  theme(panel.grid.major =element_blank(),#text=element_text(family="Times New Roman"),
        #title = element_text(family="Times New Roman"),
        #panel.grid.minor = element_blank(),
        axis.text.x= element_text(size=10, color="black"),axis.text.y= element_text(size=10),
        axis.title.x = element_text(size=12), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),strip.background = element_blank(),strip.text = element_text(size = 10))+
  geom_hline(yintercept = 1,color="#7aaed8")
ggsave(filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"ImmuOverlapInfiltrate.Odds.pdf")),
       device = "pdf",width = 5,height = 9)

x = as.data.frame(table(df$OR1>1&df$P<0.05,df$Type)[2,]/32 )
colnames(x) = "Sig"
x$no = 1-x$Sig
x$Type = row.names(x)
x %<>% pivot_longer(.,cols=1:2)
ggplot(x, aes(x = 1, y = value, fill = name)) +
  geom_bar(stat = "identity", width = 1) +    
  coord_polar(theta = "y")+facet_wrap(.~Type,ncol = 1)+theme_classic()+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size = 12),strip.text.y = element_text(size = 12,angle=0),
        strip.background = element_blank(),strip.text.y.left = element_text(angle = 0),
        legend.position = "bottom",legend.direction = "horizontal",legend.title = element_blank())+
  geom_text(aes(label = round(value*100,2)), position = position_stack(vjust = 0.5))+
  scale_fill_manual(values =c("#FDECE0","#FF6600"))
ggsave(filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"ImmuOverlapInfiltrate.Odds.Freq.pdf")),
       device = "pdf",width = 4,height = 9)

m = count[rev(1:nrow(count)),seq(1,26,5)]
colnames(m) %<>% lapply(., function(x)strsplit(x,split = ".",fixed = TRUE)[[1]][1]) %>% 
  unlist() %>% gsub("B","B cell",.) %>% gsub("CD4","CD4 cell",.)%>% gsub("CD8","CD8 cell",.)

pdf(file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"ImmuOverlapInfiltrate.pdf")),width = 5,height = 9)
corrplot(as.matrix(m), 
         "pie", 
         cl.pos = "b",
         tl.pos = "lt",
         tl.col = "black",
         is.corr = F,col=colorRampPalette(c("white", "orange"))(50),)
dev.off()
