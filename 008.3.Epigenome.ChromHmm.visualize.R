###################
#Freq.region
rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
rgb2hex(255,0,0)
rm(list = ls());gc();rm(list = ls())
Num = "008.3."
if (FALSE) { #execute
  a = fread("~/Pseudo/Result/Roadmap/Chromhmm/all.mnemonics.bed",header = FALSE)
  
  gene = fread("~/Pseudo/Result/Roadmap/Chromhmm/all.mnemonics.gene.bed",header = FALSE)
  gene %<>% dplyr::filter(.,grepl("pseudogene|protein_coding|lncRNA",V10,ignore.case = T))
  gene$V10 = gsub("protein_coding","Coding",gene$V10)
  gene[gene$V10!="Coding" & gene$V10!="lncRNA","V10"] = "Non-dynamic pseudogene"
  
  ddg = fread("~/Pseudo/Result/Human/Savedata/DDG/all.ddg.csv",header = FALSE)
  gene[gene$V9 %in% ddg$V1,"V10"] = "Dynamic pseudogene"
  gene$epitype = paste(gene$V5,gene$V10,sep = "xYx")
  typesum = tapply(gene$V13, gene$epitype, sum) %>% as.data.frame()
  #typesum$type = row.names(typesum)
  typesum$region = lapply(row.names(typesum),function(x)strsplit(x,"xYx",fixed=T)[[1]][2]) %>% unlist()
  typesum$epitype = lapply(row.names(typesum),function(x)strsplit(x,"xYx",fixed=T)[[1]][1]) %>% unlist()
  typesum = typesum[order(typesum$genetype),]
  for (i in 1:4) {
    typesum[((i-1)*15+1):(i*15),1] = typesum[((i-1)*15+1):(i*15),1] / sum(typesum[((i-1)*15+1):(i*15),1])
  }
  
  b = tapply(a$V4, a$V5, sum) %>% as.data.frame()
  b$. = b$. / sum(b$.)
  b$region = "Genome"
  b$epitype = row.names(b)
  
  b = rbind(typesum,b)
  b$num = lapply(row.names(b),function(x)strsplit(x,"_",fixed=T)[[1]][1]) %>% unlist() %>% as.numeric()
  b = b[order(b$num),]
  b$epitype %<>% factor(.,levels = unique(.))
  b$region = factor(b$region,levels = c("Genome","Coding","Dynamic pseudogene","Non-dynamic pseudogene","lncRNA"))
  b$Feature = "chromHMM"
  write.table(b,file = "~/Pseudo/Result/Roadmap/Chromhmm/Freq.region.txt",
              quote = FALSE,row.names = FALSE,sep = "\t",col.names = TRUE)
}
b = read.csv(file = "~/Pseudo/Result/Roadmap/Chromhmm/Freq.region.txt",header =TRUE,sep = "\t")
#b %<>% dplyr::filter(.,region!="lncRNA")
b = b[order(b$num),]
b$epitype %<>% factor(.,levels = unique(.))
b$region %<>% gsub("Dynamic pseudogene","Dynamic\npseudogene",.) %>% 
  gsub("Non-dynamic pseudogene","Non-dynamic\npseudogene",.)
b$region = factor(b$region,levels = c("Coding","lncRNA","Dynamic\npseudogene","Non-dynamic\npseudogene","Genome"))
p = ggplot(b,aes(region,.))+geom_bar(aes(fill=epitype),stat = "identity",position = "fill")+
  scale_fill_manual(values=c("#FF0000","#FF4500","#32CD32","#008000","#006400","#C2E105","#FFFF00","#66CDAA","#8A91D0","#CD5C5C","#E9967A","#BDB76B","#808080","#C0C0C0","#e5e5e5"))+
  theme_bw()+ylab("Proportion in state")+
  theme(axis.title.x=element_blank(),legend.position="bottom",
        axis.text.x = element_text(size=12,angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
  labs(fill="State")
ggsave(p,filename = file.path("~/Pseudo/Result/Roadmap/Chromhmm/Picture",paste0(Num,"Freq.region.pdf")),
       device = "pdf",width = 6.8,height = 6)

#####################
#Dynamic distribution
rm(list = ls());gc();rm(list = ls())
Num = "008.3."
gene = fread("~/Pseudo/Result/Roadmap/Chromhmm/all.mnemonics.gene.bed",header = FALSE) %>% as.data.frame()
n = paste(gene$V5,gene$V9,sep = ".") %>% unique() %>% as.data.frame()
n$gene = lapply(n$., function(x)strsplit(x,".",fixed=T)[[1]][2]) %>% unlist()
freq = table(n$gene) %>% as.data.frame()
freq$type = gene[match(freq$Var1,gene$V9),10]
freq = dplyr::filter(freq,grepl("pseudogene|protein_coding|lncRNA",type,type,ignore.case = T))
freq$type %<>% gsub("protein_coding","Coding",.)
freq[freq$type!="Coding" & freq$type!="lncRNA","type"] = "Non-dynamic\npseudogene"
rm(gene)
ddg = fread("~/Pseudo/Result/Human/Savedata/DDG/all.ddg.csv",header = FALSE)
freq[freq$Var1 %in% ddg$V1,"type"] = "Dynamic\npseudogene"
freq$type %<>% factor(.,levels = c("Coding","lncRNA","Dynamic\npseudogene","Non-dynamic\npseudogene"))
wilcox.test(freq[freq$type=="Dynamic\npseudogene",2],freq[freq$type=="Non-dynamic\npseudogene",2])
ggplot(freq,aes(Freq))+geom_density(aes(color=type),adjust=2)+
  xlab("Number of ChromHMM states")+ylab("Proportion")+theme_bw()+ 
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_text(size=14),legend.position='none',
    axis.text.x = element_text(size=12),
    axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
  scale_color_manual(values = c("#d95f0d","#9ecae1","#61439A","#4F69B5"))+
  coord_cartesian(xlim = c(0,15))
ggsave(filename = file.path("~/Pseudo/Result/Roadmap/Chromhmm/Picture",paste0(Num,"Proportion.number.pdf")),
       device = "pdf",width = 4.5,height = 4)
  
# ggplot(freq,aes(Freq))+geom_bar(aes(y = ..prop..,group=type,fill=type))+
#   xlab("Number of ChromHMM states")+ylab("Proportion")+theme_bw()+ 
#   theme(#panel.grid.major = element_blank(),
#         #panel.grid.minor = element_blank(),
#         axis.title.x=element_text(size=14),legend.position='none',
#         axis.text.x = element_text(size=12),
#         axis.title.y=element_text(size=14),axis.text.y = element_text(size=12))+
#   scale_fill_manual(values = c("#EE1F26","#61439A","#4F69B5","#862461"))+
#   scale_y_continuous(expand = c(0, 0),limits = c(0,.4))
# ggsave(filename = file.path("~/Pseudo/Result/Roadmap/Chromhmm/Picture",paste0(Num,"Proportion.number.pdf")),
#        device = "pdf",width = 4.5,height = 4)

## add age
Num = "008.3"
species = "Human"
age = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),header = TRUE,sep = ",",stringsAsFactors = FALSE)
freq$age = age[match(freq$Var1,age$Human),"age"]
freq %<>% na.omit(.)
cor = cor.test(freq$Freq,freq$age,method = "spearman")
freq$age %<>% as.factor() %>% as.numeric()
p = ggplot(freq,aes(age,Freq,group=age,color=age))+geom_boxplot(outlier.colour = "white",notch = TRUE)+
  geom_smooth(se=TRUE, aes(group=1),method = "lm")+
  scale_color_gradient(high = "#D6B6BB",low = "#9c0b1f")+
  theme_bw()+xlab("From old to young")+ylab("Number of state")+
  theme(axis.text.x = element_blank(),axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),axis.title.y = element_text(size=16),
        axis.ticks.x = element_blank(),panel.grid.major =element_blank(), panel.grid.minor = element_blank())+
  guides(color=FALSE)+
  annotate(geom="text", x=8.5, y=13.5, size=6,
           label=paste0("rho=",round(cor$estimate,2),",","\n","P=",round(cor$p.value,4)))
ggsave(p, filename = paste0("~/Pseudo/Result/Roadmap/Chromhmm/Picture/",Num,".state.age.pdf"),
       width = 5,height = 4)

## add gene type
S=species
b = read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
         header = FALSE,sep = "\t",stringsAsFactors = F)
m1 <- dplyr::filter(b,grepl("polymorphic_pseudogene",V5))
m2 <- dplyr::filter(b,grepl("unprocessed_pseudogene",V5))
m3 <- dplyr::filter(b,grepl("unitary_pseudogene",V5))
m4 <- dplyr::filter(b,V4 %in% setdiff(b$V4,c(m1$Gene,m2$Gene,m3$Gene)))
m <- rbind(m1,m2,m3,m4)
m$type2 = rep(c("Polymorphic","Unprocessed","Unitary","Processed"),c(nrow(m1),nrow(m2),nrow(m3),nrow(m4)))
freq$type = m[match(freq$Var1,m$V4),"type2"]

options(digits = 2)
df = matrix(0,nrow = 4,ncol = 2) %>% as.data.frame()
for (i in 1:4) {
  cor = with(freq[freq$type==unique(freq$type)[i],],cor.test(Freq,as.numeric(age),method = "spearman"))
  df[i,1] = cor$estimate
  df[i,2] = cor$p.value
}
row.names(df) = unique(freq$type)
df$V3 = paste0("rho=",round(df$V1,2),",","\n","P=",format(df$V2,2))

#with(freq[freq$type=="Poly",],cor.test(Freq,age,method = "spearman"))
freq$age %<>% as.factor()
p = ggplot(freq,aes(age,Freq))+geom_boxplot(outlier.colour = "white",notch = TRUE)+
  geom_smooth(se=TRUE, aes(group=1),method = "lm")+
  scale_color_gradient(high = "#D6B6BB",low = "#9c0b1f")+
  theme_bw()+xlab("From old to young")+ylab("Number of state")+
  theme(axis.text.x = element_blank(),axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=14),axis.title.y = element_text(size=16),
        axis.ticks.x = element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12),strip.background = element_rect(fill="white"))+
  guides(color=FALSE)+
  #annotate(geom="text", x=8.5, y=13.5, size=6,
  #         label=paste0("rho=",round(cor$estimate,2),",","\n","P=",round(cor$p.value,4)))+
  facet_wrap(.~type,nrow = 2)+
  geom_text(aes(x,y,label =lab),
            data = data.frame(x = 8.5,
                              y = 15,
                              lab = as.character(df$V3),
                              #sex = factor(c("brain","liver","testis")),levels = c("brain","liver","testis")),
                              type = factor(unique(freq$type)),levels = unique(freq$type)), #facet use name
            vjust = 1, size=4)
ggsave(p, filename = paste0("~/Pseudo/Result/Roadmap/Chromhmm/Picture/",Num,".state.age.type.pdf"),
       width = 7,height = 6)
