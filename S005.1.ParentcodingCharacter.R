rm(list = ls());gc();rm(list = ls())
Num = "S005.1."
#cd ~/Pseudo/Data/Seqdata/EXAC
#wget https://storage.googleapis.com/gnomad-public/legacy/exacv1_downloads/release0.3.1/manuscript_data/forweb_cleaned_exac_r03_march16_z_data_pLI.txt.gz
a = read.csv("~/Pseudo/Data/Ref/Human/pseudo.coding.pair.txt",header = T,sep = "\t")
S = "Human"
b = read.csv(file.path("/home/qians/Pseudo/Data/Ref",S,paste0("gene",S,".bed")),
             header = FALSE,sep = "\t",stringsAsFactors = F)
b %<>% dplyr::filter(.,V1 %in% c(1:100,"X"))
a$type = b[match(a$Pseudogene,b$V4),5]
pli = fread("~/Pseudo/Data/Seqdata/EXAC/forweb_cleaned_exac_r03_march16_z_data_pLI.txt") %>% as.data.frame()
id = fread("~/Pseudo/Data/Ref/Human/gene.ENSG.name.2.txt",header = FALSE) %>% as.data.frame()

# pLI ~ Number ----------------
#* Processed -----------------------------------------------------------
tmp = dplyr::filter(a,grepl("processed",type)) %>% dplyr::filter(.,!grepl("unprocessed",type))
freq = table(tmp$Parentgene) %>% as.data.frame()

freq$id = id[match(freq$Var1,id$V1),2]
freq$pli = pli[match(freq$id,pli$gene),"pLI"]
freq %<>% na.omit()
cor = cor.test(freq$Freq,freq$pli,method = "spearman")
ggplot(freq,aes(pli,log10(Freq)))+geom_point()+geom_smooth(method = "lm")+
  theme_bw()+xlab("Probability of being LoF intolerant (pLI)")+
  ylab("Number of processed pseudogenes")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=12),axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=12),axis.title.y = element_text(size=14))+
  annotate(geom="text", x=0.15, y=2.05, size=5,
           label=paste0("rho=",round(cor$estimate,2),",","\n","P=",format(cor$p.value,digits = 2)))
#freq$group = cut(freq$pli,breaks = seq(0,1,0.1))

#* Unprocessed ---------------------------------------------------------
tmp = dplyr::filter(a,grepl("unprocessed",type)) 
freq = table(tmp$Parentgene) %>% as.data.frame()

freq$id = id[match(freq$Var1,id$V1),2]
freq$pli = pli[match(freq$id,pli$gene),"pLI"]
freq %<>% na.omit()
cor.test(freq$Freq,freq$pli,method = "spearman")
ggplot(freq,aes(pli,log10(Freq)))+geom_point()+geom_smooth(method = "lm")

#* Expressed ---------------------------------------------------------------
exp = fread("~/Pseudo/Result/Human/Savedata/gene.expressnum.0.3allfpkm.csv") %>% as.data.frame()
a$exp = exp[match(a$Pseudogene,exp$gene),2]
tmp = dplyr::filter(a, exp > 1) 
freq = table(tmp$Parentgene) %>% as.data.frame()

freq$id = id[match(freq$Var1,id$V1),2]
freq$pli = pli[match(freq$id,pli$gene),"pLI"]
freq %<>% na.omit()
cor.test(freq$Freq,freq$pli,method = "spearman")
ggplot(freq,aes(pli,log10(Freq)))+geom_point()+geom_smooth(method = "lm")


#* Chr ---------------------------------------------------------------------
freq = table(a$Parentgene) %>% as.data.frame()
freq$chr = b[match(freq$Var1,b$V4),1]
freq %<>% dplyr::filter(.,chr %in% c(1:100,"X"))
table(b[grepl("coding",b$V5),1]=="X")/tapply(freq$Freq, freq$chr=="X", sum)

#* Chr + Processed ---------------------------------------------------------------------
tmp = dplyr::filter(a,grepl("processed",type)) %>% dplyr::filter(.,!grepl("unprocessed",type))
freq = table(tmp$Parentgene) %>% as.data.frame()
freq$chr = b[match(freq$Var1,b$V4),1]
freq %<>% dplyr::filter(.,chr %in% c(1:100,"X"))
wil = wilcox.test(freq[freq$chr=="X",2],freq[freq$chr!="X",2])
p = ggplot(freq,aes(log(Freq),color=chr))+geom_density()+
  theme_bw()+xlab("Number of processed pseudogene (log10)")+ylab("Frequency")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=12),axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=12),axis.title.y = element_text(size=14))+
  scale_color_manual(values = c(rep("#8db8d5",length(unique(freq$chr))-1),"#9c0b1f"))+
  geom_vline(xintercept = median(freq[freq$chr!="X","Freq"]),color="#035782",linetype="dashed")+
  geom_vline(xintercept = median(freq[freq$chr=="X","Freq"]),color="#9c0b1f",linetype="dashed")+
  guides(color=FALSE)+
  annotate("text",x=1.5,y=0.5,label=paste0("P=",format(wil$p.value,digits = 2)),size=5)
ggsave(p, filename = file.path("/home/qians/Pseudo/Result/Human/Picture",paste0(Num,"Frequency.Processed.AX.pdf")),
       width = 5,height = 4)

#* C+E+P ---------------------------------------------------------------------
exp = fread("~/Pseudo/Result/Human/Savedata/gene.expressnum.0.3allfpkm.csv") %>% as.data.frame()
a$exp = exp[match(a$Pseudogene,exp$gene),2]
tmp = dplyr::filter(a,grepl("processed",type)) %>% 
  dplyr::filter(.,!grepl("unprocessed",type)) %>% 
  dplyr::filter(., exp > 1) 
freq = table(tmp$Parentgene) %>% as.data.frame()
freq$chr = b[match(freq$Var1,b$V4),1]
freq %<>% dplyr::filter(.,chr %in% c(1:100,"X"))
wil = wilcox.test(freq[freq$chr=="X",2],freq[freq$chr!="X",2])
p = ggplot(freq,aes(log(Freq),color=chr))+geom_density()+
  theme_bw()+xlab("Number of expressed processed pseudogene (log10)")+ylab("Frequency")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=12),axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=12),axis.title.y = element_text(size=14))+
  scale_color_manual(values = c(rep("#8db8d5",length(unique(freq$chr))-1),"#9c0b1f"))+
  geom_vline(xintercept = median(freq[freq$chr!="X","Freq"]),color="#035782",linetype="dashed")+
  geom_vline(xintercept = median(freq[freq$chr=="X","Freq"]),color="#9c0b1f",linetype="dashed")+
  guides(color=FALSE)+
  annotate("text",x=1.5,y=0.5,label=paste0("P=",format(wil$p.value,digits = 2)),size=5)
ggsave(p,filename = file.path("/home/qians/Pseudo/Result/Human/Picture",paste0(Num,"Frequency.ExpressedProcessed.AX.pdf")),
       width = 5,height = 4)

# Express ~ Number ----------------
## ALL
exp = fread("~/Pseudo/Data/Seqdata/Illumina/Human/Human.allgene.fpkm.txt") %>% 
  as.data.frame()
exp$type1 = "No-Parent"
exp[exp$V1 %in% a$Parentgene,"type1"] = "Parent"
exp$type2 = b[match(exp$V1,b$V4),5]
exp %<>% na.omit() %>% dplyr::filter(.,grepl("coding",type2))
df = t(apply(exp[,2:314], 2, function(x)tapply(x,exp$type1,median))) %>% as.data.frame()
df$ratio = df$Parent / df$`No-Parent`
df$type ="A"
ggplot(df,aes(type,ratio))+geom_boxplot(notch = TRUE,outlier.colour = "white",fill="#00434f")+
  theme_classic()+ylab("Expression ratio")+xlab("Parent/No-parent")+
  #theme_bw()
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=12),axis.title.y = element_text(size=14))+
  coord_cartesian(ylim = c(0,7))+
  annotate("pointrange", x = 1, y = 1, ymin = 1, ymax = 1,
           colour = "#ec0b05", size = 0.5)+
  annotate("text", x=1, y=6.5, label="***",size=6)
ggsave(filename = file.path("/home/qians/Pseudo/Result/Human/Picture",paste0(Num,"ExpressionRatio.ParentNoParent.pdf")),
       width = 2.5,height = 3)
  
## Processed
exp = fread("~/Pseudo/Data/Seqdata/Illumina/Human/Human.allgene.fpkm.txt") %>% 
  as.data.frame()
exp$type1 = "No-Parent"
exp[exp$V1 %in% a[grepl("process",gsub("unprocessed","AAA",a$type)),"Parentgene"],"type1"] = "Parent"
exp$type2 = b[match(exp$V1,b$V4),5]
exp %<>% na.omit() %>% dplyr::filter(.,grepl("coding",type2))
df = t(apply(exp[,2:314], 2, function(x)tapply(x,exp$type1,median))) %>% as.data.frame()
df$ratio = df$Parent / df$`No-Parent`
df$type ="A"
ggplot(df,aes(type,ratio))+geom_boxplot(notch = TRUE,outlier.colour = "white",fill="#00434f")+
  theme_classic()+ylab("Expression ratio")+xlab("Parent/No-parent")+
  #theme_bw()
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=12),axis.title.y = element_text(size=14))+
  coord_cartesian(ylim = c(0,8))+
  annotate("pointrange", x = 1, y = 1, ymin = 1, ymax = 1,
           colour = "#ec0b05", size = 0.5)+
  annotate("text", x=1, y=7, label="***",size=6)
ggsave(filename = file.path("/home/qians/Pseudo/Result/Human/Picture",paste0(Num,"ExpressionRatio.ParentNoParent.processed.pdf")),
       width = 2.5,height = 3)


