rm(list = ls());gc();rm(list = ls())
# GUSBP1 in dissertation  (search QuaDMutEx gusbp1)
# https://scholarscompass.vcu.edu/cgi/viewcontent.cgi?article=6724&context=etd
Num = "S003.1."
a = fread("~/Pseudo/Data/Seqdata/CaiXue2020natureRICseq/HeLa_merge.intermolecular_interaction.used.bed") %>%
  as.data.frame()

#b = a[a$geneA_symbol=="SUZ12P"|a$geneB_symbol=="SUZ12P",]

x = a[a$geneA_symbol=="SUZ12P"|a$geneB_symbol=="SUZ12P",]
x = x[x$geneA_type=="protein_coding"|x$geneB_type=="protein_coding",]
x = c(x$geneA_symbol,x$geneB_symbol) %>% unique()

#Number of interaction
if (TRUE) {
  freq1 = table(a$geneA_symbol) %>% as.data.frame()
  #freq1$type = a[match(freq1$Var1,a$geneA_symbol),"geneA_type"]
  freq2 = table(a$geneB_symbol) %>% as.data.frame()
  #freq2$type = a[match(freq2$Var1,a$geneB_symbol),"geneB_type"]
  
  freq = rbind(freq1,freq2)
  freq = tapply(freq$Freq, freq$Var1, sum) %>% as.data.frame()
  
  type = fread("~/Pseudo/Data/Ref/Human/hg19.gene.ens.id.type.txt",header = FALSE) %>% as.data.frame()
  freq[,c("ens","type")] = type[match(row.names(freq),type$V2),c(1,3)]
  freq %<>% na.omit() %>%dplyr::filter(.,grepl("protein|pseudo|lincRNA",type))
  table(freq$type)
  #freq$type2 = gsub("transcribed_","",freq$type) %>% gsub("_pseudogene","",.) %>% 
  #  gsub("translated_","",.) %>% gsub("protein_","",.)
  freq$type2 = gsub("IG_V_","",freq$type) %>% gsub("TR_V_","",.) %>% 
    gsub("polymorphic_","",.) %>% gsub("protein_","",.)
  
  #freq %<>% dplyr::filter(.,type2 %in% c("coding","unprocessed","processed","unitary","polymorphic"))
  
  ddg = fread("~/Pseudo/Result/Human/Savedata/DDG/all.ddg.csv",header = FALSE) %>% as.data.frame()
  freq$type3 = freq$type2 %>% gsub("coding","Coding",.) %>% gsub("lincRNA","lncRNA",.) %>% gsub("pseudogene","Non-dynamic",.)
  #freq[freq$type2!="coding","type3"] = "Nondynamic"
  freq[freq$ens %in% ddg$V1,"type3"] = "Dynamic"
  table(freq$type3)
  #plot
  ratio = table(freq$type3) %>% as.data.frame()
  ratio$all = c(20327,sum(unique(ddg$V1) %in% type$V1),7109,
                nrow(type[grepl("pseu",type$V3),])-sum(unique(ddg$V1) %in% type$V1))
  ratio$ratio = ratio$Freq / ratio$all
  
  ratio[c(1,3),2:3] %>% fisher.test()
  ratio[2:3,2:3] %>% fisher.test()
  ratio[c(3,4),2:3] %>% fisher.test()
  ratio$Var1 %<>% factor(.,levels = c("Coding","lncRNA","Dynamic","Non-dynamic"))
  ggplot(ratio,aes(Var1,ratio,fill=Var1))+geom_bar(stat = "identity",position = "dodge")+
    theme_classic()+ylab("Ratio of trans interacted")+
    theme(axis.text.x= element_text(size=12),axis.text.y= element_text(size=12),legend.position='none',
          axis.title.x = element_blank(), axis.title.y = element_text(size=14))+
    scale_fill_manual(values =  c("#d95f0d","#9ecae1","#61439A","#4F69B5"))+
    scale_y_continuous(expand = c(0, 0),limits = c(0,.9))+
    geom_segment(aes(x=1, y=0.85, xend=2, yend=0.85))+
    annotate("text", x=1.5, y=0.86, label="***",size=5)+
    geom_segment(aes(x=2, y=0.48, xend=3, yend=0.48))+
    annotate("text", x=2.5, y=0.51, label="n.s.",size=5)+
    geom_segment(aes(x=3, y=0.45, xend=4, yend=0.45))+
    annotate("text", x=3.5, y=0.46, label="***",size=5)
  ggsave(filename = file.path("~/Pseudo/Result/Human/Picture",paste0(Num,"Ratio.transInteracted.pdf")),
         device = "pdf",width = 4.5,height = 3.7)
  
  my_comparisons <- list(c("Coding","Dynamic"),c("Dynamic","Nondynamic"))
  ggplot(freq,aes(type3,log10(.),fill=type3))+geom_boxplot(notch = TRUE,outlier.colour = "white")+
    theme_classic()+ylab("Num of trans interacted (log10)")+
    theme(axis.text.x= element_text(size=12),axis.text.y= element_text(size=12),legend.position='none',
          axis.title.x = element_blank(), axis.title.y = element_text(size=14))+
    coord_cartesian(ylim = c(0,3))+
    scale_fill_manual(values = c("#EE1F26","#61439A","#4F69B5"))+
    stat_compare_means(comparisons = my_comparisons,label.y = c(2.7,2.5),tip.length = 0.015)
  ggsave(filename = file.path("~/Pseudo/Result/Human/Picture",paste0(Num,"Num.transInteracted.pdf")),
         device = "pdf",width = 3.6,height = 3.7)

b = a[,c(4,5,10,11,13,15)]
id = fread("~/Pseudo/Data/Ref/Human/hg19.gene.ens.id.type.txt",header = FALSE) %>% as.data.frame()
b$ida = id[match(b$geneA_symbol,id$V2),1]
b$idb = id[match(b$geneB_symbol,id$V2),1]

dim(a[a$geneA_symbol=="SUZ12P"|a$geneB_symbol=="SUZ12P",])


b %<>% na.omit() %>% dplyr::filter(.,grepl("coding|pseu",geneA_type)) %>%
  dplyr::filter(.,grepl("coding|pseu",geneB_type))

ddg = fread("~/Pseudo/Result/Human/Savedata/DDG/all.ddg.csv",header = FALSE) %>% as.data.frame()
b$typea = "C"
b[b$geneA_type!="protein_coding","typea"] = "P"
b[b$ida %in% ddg$V1,"typea"]= "D"
b$typeb = "C"
b[b$geneB_type!="protein_coding","typeb"] = "P"
b[b$idb %in% ddg$V1,"typeb"]= "D"
b$type = paste0(b$typea,b$typeb)  %>% gsub("DC","CD",.) %>% gsub("PC","CP",.) %>% gsub("PD","DP",.)
ggplot(b,aes(type,log10(NumberOfChimericReads)))+geom_boxplot()

ggplot(b,aes(log(NumberOfChimericReads),-log10(Pvalue)))+geom_point()+geom_smooth()

