rm(list = ls())
Num = "001.3."

# Chr distribution (age) ----
if (TRUE) {
  rm(list = ls())
  mtx = c()
  Num = "001.3."
  for (species in c("Human","Mouse")) {
    a = read.csv(paste0("/home/qians/Pseudo/Result/",species,"/Savedata/",species,".gene.age.csv"),header = TRUE,sep = ",",stringsAsFactors = FALSE)
    a[a$chr!="X","chr"] = "A"
    freq = table(a$Max,a$chr)[,2]/c(table(a$Max,a$chr)[,1]+table(a$Max,a$chr)[,2]) %>% as.data.frame()
    freq$Var1 = 1:nrow(freq)
    freq$age = a[match(freq$Var1,a$Max),"age"]
    freq$species = species
    mtx = rbind(mtx,freq)
  }
  mtx$. = mtx$. * 100
   p = ggplot(mtx,aes(age,.))+geom_point(aes(shape=species),size=2)+geom_line(aes(linetype=species))+theme_bw()+
    xlab("Million year")+ylab("Proportion of X-linked genes (%)")+
    theme(legend.position = c(0.2, 0.9),legend.title = element_blank(),
          legend.background = element_blank(),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size=10),axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=10),axis.title.y = element_text(size=12))+
    #guides(linetype=guide_legend(title=NULL))+theme(legend.position=c(.15,.85))+
    scale_shape_manual(values = c(19,1))
  ggsave(p, filename = paste0("/home/qians/Pseudo/Result/Human/Picture/",Num,"HumanMouse.proportionX.age.pdf"),
         width = 3.3,height = 3.3)
}