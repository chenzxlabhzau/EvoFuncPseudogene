rm(list = ls());gc();rm(list = ls())
Num = "015.4."

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

count = matrix(0,nrow = length(all),ncol = 2) %>% as.data.frame()
colnames(count) = c("Type","Number")
count$Type = all
for (i in 1:length(all)) {
  focal.cancer = paste0("TCGA-",all[i])
  a = fread(file =  paste0("/home/qians/Pseudo/Data/Seqdata/TCGA/ImmunoPseudo/",
                           focal.cancer,".ImmuPseudo.csv")) %>% 
    as.data.frame() %>% na.omit()
  
  count[i,2] = sum(a$padj<0.05 & abs(a$sigValue)>0.995)
}
ggplot(count,aes(Type,Number))+geom_bar(stat = "identity",position = "dodge",fill="#FF6600")+
  coord_flip()+theme_bw()+
  ylab("Immune-related pseudogene")+
  theme(panel.grid.major =element_blank(),#text=element_text(family="Times New Roman"),
        #title = element_text(family="Times New Roman"),
        #panel.grid.minor = element_blank(),
        axis.text.x= element_text(size=12, color="black"),axis.text.y= element_text(size=12),
        axis.title.x = element_text(size=14), axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+
  scale_y_continuous(expand = c(0, 0),limits = c(0,5300))
ggsave(filename = file.path("~/Pseudo/Result/TCGA/Picture",paste0(Num,"Number.Immunerelated.pdf")),
       device = "pdf",width = 3.2,height = 5)

  
  