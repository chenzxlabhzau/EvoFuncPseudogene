rm(list = ls());gc();rm(list = ls())
Num = "006.6"

species = "Human"
df = c()
text = c()
for (species in c("Human","Mouse")) {
  a = read.csv(paste0("~/Pseudo/Result/",species,"/Savedata/DDG/all.ddg.csv"),header = FALSE)
  freq = table(a$V1) %>% as.data.frame()
  a = unique(a$V1) %>% as.data.frame()
  
  age = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),header = TRUE,sep = ",",stringsAsFactors = FALSE)
  age$num = freq[match(age[,1],freq$Var1),2]
  age[is.na(age$num),"num"] = 0
  cor.test(table(age$num>0,age$age)[2,]/table(age$age),1:10,method = "spearman")
  c = as.data.frame(table(age$num>0,age$age)[2,]/table(age$age)) 
  ggplot(c,aes(Var1,Freq))+geom_bar(stat = "identity",position = "dodge")+geom_smooth(aes(group=1),se = FALSE,method = "lm")
  c$Var1 %<>% as.character() %>% as.numeric()
  c$Freq = c$Freq * 100
  c$species = species
  df = rbind(df,c)
  cor = cor.test(c$Var1,c$Freq,method = "spearman")
  text = c(text, paste0("rho= ",round(cor$estimate,2),",","\n","P= ",round(cor$p.value,4)))
}
p = ggplot(df,aes(Var1,Freq))+geom_point()+geom_smooth(aes(group=1),se = FALSE,method = "lm")+
  theme_bw()+xlab("Million year")+ylab("Proportion of dynamic expression (%)")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.title.x=element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.title.y=element_text(size=12),
        axis.text.y = element_text(size=10),
        strip.text = element_text(size=10),strip.background = element_rect(fill="white"))+
  geom_text(aes(x,y,label =lab),
            data = data.frame(x = -80,
                              y = 60,
                              lab = text,
                              species = factor(c("Human","Mouse"),levels = c("Human","Mouse"))),
            vjust = 1, size=4)+
  facet_grid(.~species)
ggsave(p,filename = paste0("/home/qians/Pseudo/Result/Human/Picture/",Num,".HumanMouse.proportiondynamic.age.pdf"),
       width = 6,height = 3.5)
