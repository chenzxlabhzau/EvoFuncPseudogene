rm(list = ls())
Num = "002.7."

mtx = c()
text = c()
for (species in c("Human","Mouse")) {
express = fread(file.path("~/Pseudo/Result",species,"Savedata","gene.expressnum.1allfpkm.csv")) %>% as.data.frame()
age = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),header = TRUE,sep = ",",stringsAsFactors = FALSE)
express[,c("age","chr")] = age[match(express$gene,age[,1]),c("age","chr")]
express %<>% na.omit()
express[express$chr!="X","chr"]="A"
df = as.data.frame(tapply(express$num >0,express$age, sum)/table(express$age),stringsAsFactors=FALSE)
df$Var1 %<>% as.numeric()
df$species = species
cor = cor.test(df$Var1,df$Freq,method = "spearman")
text = c(text, paste0("rho= ",round(cor$estimate,2),",","\n","P= ",round(cor$p.value,4)))
mtx = rbind(mtx,df)
}
p = ggplot(mtx,aes(Var1,Freq*100))+geom_point()+facet_grid(.~species)+geom_smooth(method = "lm",se = FALSE)+
  xlab("Millon year")+ ylab("Fraction of expressed pseudogenes")+theme_bw()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.title.x=element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.title.y=element_text(size=12),
        axis.text.y = element_text(size=10),
        strip.text = element_text(size=12),
        strip.background = element_rect(fill="white"))+
  geom_text(aes(x,y,label =lab),
            data = data.frame(x = -380,
                              y = 63,
                              lab = text,
                              species = factor(c("Human","Mouse"),levels = c("Human","Mouse"))),
            vjust = 1, size=5)
ggsave(p,filename = paste0("/home/qians/Pseudo/Result/Human/Picture/",Num,"HumanMouse.proportionexpressed.age.pdf"),
       width = 7,height = 4)

p2 = mtx %>% dplyr::filter(species== "Human") %>% ggplot(aes(Var1,Freq*100))+geom_point()+geom_smooth(method = "lm",se = FALSE)+
  xlab("Millon year")+ ylab("Ratio of expressed pseudogenes (%)")+theme_classic()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.title.x=element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.title.y=element_text(size=12),
        axis.text.y = element_text(size=10),
        strip.text = element_text(size=12),
        strip.background = element_rect(fill="white",color="white"))+
  geom_text(aes(x,y,label =lab),
            data = data.frame(x = -380,
                              y = 25,
                              lab = text[1],
                              species = factor(c("Human"),levels = c("Human"))),
            vjust = 1, size=4)
ggsave(p2,filename = paste0("/home/qians/Pseudo/Result/Human/Picture/",Num,"Human.proportionexpressed.age.pdf"),
       width = 4.5,height = 3.5)

if (FALSE) {
  ggplot(express,aes(age,log(num+1)))+geom_boxplot(aes(group=age))
  ggplot(df,aes(Var1,Freq))+geom_point()
  
  #different group: A X 
  df = as.data.frame(tapply(express$num >0,paste(express$age,express$chr,sep = "_"), sum)/table(paste(express$age,express$chr,sep = "_")))
  df$age = df$Var1 %>% lapply(.,function(x)strsplit(as.character(x),split = "_",fixed = TRUE)[[1]][1]) %>% unlist() %>% as.numeric()
  df$chr = df$Var1 %>% lapply(.,function(x)strsplit(as.character(x),split = "_",fixed = TRUE)[[1]][2]) %>% unlist() 
  ggplot(df,aes(age,Freq))+geom_point()+facet_wrap(.~chr)+geom_smooth()
}
