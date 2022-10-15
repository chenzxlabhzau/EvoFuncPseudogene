{
  rm(list = ls())
  S= c("Human","Mouse")
  dir1 = file.path("/home/qians/Pseudo/Result",S[1],"Savedata/DDG")
  file1 = grep("ddg.csv$",list.files(dir1),value=TRUE)
  dir2 = file.path("/home/qians/Pseudo/Result",S[2],"Savedata/DDG")
  file2 = grep("ddg.csv$",list.files(dir2),value=TRUE)
  if (all(file1==file2)) {
    j= lapply(file1, function(x)strsplit(x,"_",fixed=T)[[1]][1]) %>% unlist()
    c = matrix(0,nrow = length(S) * length(file1),ncol = 3) %>% as.data.frame()
    colnames(c) = c("Species","Tissue","Number")
    c$Species = rep(S,each= length(file1))
    c$Tissue = rep(j,2)
    c %<>% dplyr::filter(Tissue %in% c("Brain","Cerebellum","Heart","Kidney","Liver","Ovary","Testis"))
    for (i in 1:nrow(c)) {
      gene = read.csv(file.path("/home/qians/Pseudo/Result",c[i,1],"Savedata/DDG",paste0(c[i,2],"_ddg.csv")))
      c[i,3] = nrow(gene)
    }
  }
}
ggplot(c,aes(Tissue,Number))+geom_bar(aes(fill=Species),stat = "identity",position ="dodge")+
  theme_classic()+facet_grid(.~Species)+ylab("Count")+
  theme(axis.title.x=element_blank(),legend.position='none',
        axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1),
        axis.title.y=element_text(size=14),axis.text.y = element_text(size=12),
        strip.background  = element_blank(),
        strip.text = element_text(size=8))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,3000))   #+scale_fill_manual(values=c("#e8ba31","#56a4d8"))
ggsave(filename = file.path("~/Pseudo/Result","Human","Picture","DDG.pdf"),
       device = "pdf",width = 5.3,height = 3)
