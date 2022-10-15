#4wpcref

#cd ~/Pseudo/Result/Human/Savedata/DEG/4wpcref
#wc -l *txt | grep txt | awk -F " " '{print $1"\t"$2}' > 4wpcref.DEG.number

rm(list = ls())
a = read.csv("/home/qians/Pseudo/Result/Human/Savedata/DEG/4wpcref/4wpcref.DEG.number",
             header = F,sep = "\t",stringsAsFactors = F)
wd = "/home/qians/Pseudo/Data/Seqdata/Illumina/"
Species = "Human"

for (i in 1:6) {
  a[,i+2] = lapply(a$V2,function(x)strsplit(x,".",fixed=T)[[1]][i]) %>% unlist()
}
deve = read.csv(paste0(wd,Species,"/",Species,".developstage.txt"),
                header = F,sep = "\t",stringsAsFactors = F,row.names = 1)
deve["Number",] = colnames(deve) %>% gsub("V","",.) 
a$V9 = a$V5
for (i in 1:7) {
  for (j in 2:12) {
    a$V9 %<>% gsub(deve[i,j],deve[8,j],.) #
  }
}
a$V9 %<>% as.numeric()
a$V9 = a$V9-2
tapply(a$V1, a$V8, mean)

for (i in 1:nrow(a)) {
  if (a[i,8]=="coding") {
    a[i,1] = a[i,1]/10 #log10(a[i,1]+1)
  }
}

for (i in 1:nrow(a)) {
  if (a[i,7]=="down") {
    a[i,1] = -a[i,1]
  }
}

#ggplot(a,aes(V9,V1,color=V8))+geom_point()+geom_line(aes(group=V8))+facet_grid(.~V4)

#dplyr::filter(a,V4=="Brain") %>% ggplot(.,aes(V9,V1,color=V8,linetype=V7))+geom_line(aes(group=V8))
a$V7 %<>% factor(.,levels = c("up","down"))
colnames(a)[7:8] = c("DE","Type")
p = a %>% ggplot(.,aes(V9,V1,color=Type,linetype=DE))+geom_line()+
  geom_vline(xintercept = 8,linetype=2,color="grey",size=1)+
  xlab("Developmenttal")+ylab("DE genes")+
  #geom_hline(yintercept = 0)+
  theme_classic()+facet_grid(.~V4)+
  theme(legend.position = 'top', legend.direction = 'horizontal',
        axis.title.x=element_blank(),#legend.position='none',
        axis.text.x = element_text(size=18,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=20),axis.text.y = element_text(size=18),
        strip.background  = element_blank(),
        strip.text = element_text(size=12))+
  scale_x_continuous(breaks = c(3,6,9),labels = c("8wpc","13wpc","Toddler"))
ggsave(p,filename = file.path("~/Pseudo/Result",Species,"Picture",paste0(Species,".number.DEG.4wpcref.pdf")),
       device = "pdf",width = 11,height = 6)
