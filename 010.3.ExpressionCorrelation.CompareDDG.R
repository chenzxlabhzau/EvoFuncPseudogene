rm(list = ls());gc();rm(list = ls())
#cd ~/Pseudo/Data/Ref/Human
#zcat Homo_sapiens.GRCh37.87.chr.gtf.gz | grep -w gene | awk -F ";" '{print $1"\t"$3"\t"$5}' | awk -F "\t" '{print $9"\t"$10"\t"$11}'|sed 's/gene_id//g' | sed 's/gene_name//g' | sed 's/gene_biotype//g' | sed 's/"//g' | sed 's/ //g' > hg19.gene.ens.id.type.txt

##################### co-expression number
#S="Human"
wd = "/home/qians/Pseudo/Data/Seqdata/Illumina"
library(simplifyEnrichment)
library(Hmisc)

if (FALSE) { #has executed
  for (S in c("Mouse","Human")) {
    pair = fread(file = file.path("~/Pseudo/Result",S,"Savedata","expressionPair.csv"),nThread = 5)
    colnames(pair)[1] = "pseu"
    pair.filter = dplyr::filter(pair,abs(R)>=0.9 & Pvalue < 0.01) %>% as.data.frame()
    fwrite(pair.filter,file = file.path("~/Pseudo/Result",S,"Savedata","expressionPair.sig.csv"),nThread = 5)
  }
}  

for (S in c("Mouse","Human")) {  
  pair.filter = read.csv(file = file.path("~/Pseudo/Result",S,"Savedata","expressionPair.sig.csv"),header = TRUE)
  ddpseu = fread(paste0("~/Pseudo/Result/",S,"/Savedata/DDG/all.ddg.csv"),header = FALSE) 
  n = table(ddpseu$V1) %>% as.data.frame()
  pair.filter$ddpseu = 0
  pair.filter[pair.filter$pseu %in% ddpseu$V1,"ddpseu"]=1
  pair.filter$freq = n[match(pair.filter$pseu,n$Var1),2]
  pair.filter[is.na(pair.filter)]=0
  
  ddcoding = fread(file.path(wd,S,paste0(S,".gene.development.annotation.csv"))) %>% 
    as.data.frame()
  ddcoding %<>% dplyr::select(.,c(ends_with("ID"),ends_with("DDG"))) #%>% 
  # dplyr::select(.,c(ends_with("ID"),starts_with(Tissue)))
  ddcoding[is.na(ddcoding)]=0
  ddcoding$number = apply(ddcoding[,-1],1,sum)
  pair.filter$ddcoding = ddcoding[match(pair.filter$coding,ddcoding[,1]),"number"]
  
  # Enrichment of dynamic pseu
  c(sum(n$Var1 %in% pair.filter$pseu),nrow(n),length(unique(pair.filter$pseu)),
    nrow(read.csv(file.path(wd,S,paste0(S,".pseudogene.tpm.txt")),header = T,sep = "\t"))) %>% 
    matrix(.,nrow=2) %>% fisher.test()
  sum(n$Var1 %in% pair.filter$pseu)/ nrow(n); 
  (length(unique(pair.filter$pseu))-sum(n$Var1 %in% pair.filter$pseu)) / 
    (nrow(read.csv(file.path(wd,S,paste0(S,".pseudogene.tpm.txt")),header = T,sep = "\t"))-nrow(n))
  # Dynamic psuedogenes have more interactions
  inter.number = table(pair.filter$pseu) %>% as.data.frame()
  inter.number$type = pair.filter[match(inter.number$Var1,pair.filter$pseu),"ddpseu"]
  inter.number$type %<>% gsub(1,"Dynamic",.) %>% gsub(0,"Non-dynamic",.)
  wilcox.test(inter.number[inter.number$type=="Dynamic",2],
              inter.number[inter.number$type=="Non-dynamic",2],alternative = "greater")
  r=(length(unique(pair.filter$pseu))-sum(n$Var1 %in% pair.filter$pseu)) / 
    (nrow(read.csv(file.path(wd,S,paste0(S,".pseudogene.tpm.txt")),header = T,sep = "\t"))-nrow(n))
  if (S=="Mouse") {
    label = 0.23
    p = ggplot(inter.number,aes(type,log10(Freq),fill=type))+geom_boxplot(notch = T)+theme_classic()+
      ylab("Co-expressed pair number (log10)")+
      theme(axis.title.x=element_blank(), axis.title.y=element_text(size=16),
            axis.text.x = element_text(size=14,angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size=14))+
      scale_fill_manual(values =c("#61439A","#4F69B5")) +
      theme(legend.position='none')+
      geom_segment(aes(x=1, y=3, xend=2, yend=3))+
      annotate("text", x=1.5, y=3.1, label="***",size=5)+
      annotate("text", x=0.90, y=label, label=
                 sum(n$Var1 %in% pair.filter$pseu))+
      annotate("text", x=1.16, y=label, label=
                 paste0("(",round(sum(n$Var1 %in% pair.filter$pseu)/nrow(n)*100,2),"%)"))+
      annotate("text", x=1.93, y=label, label=
                 length(unique(pair.filter$pseu))-sum(n$Var1 %in% pair.filter$pseu))+
      annotate("text", x=2.16, y=label, label=
                 paste0("(",round(r,2),"%)"))
  }else{
    label = 0.25 #0.3
    p = ggplot(inter.number,aes(type,log10(Freq),fill=type))+geom_boxplot(notch = T)+theme_classic()+
      ylab("Co-expressed pair number (log10)")+
      theme(axis.title.x=element_blank(), axis.title.y=element_text(size=16),
            axis.text.x = element_text(size=14,angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size=14))+
      scale_fill_manual(values =c("#61439A","#4F69B5")) +
      theme(legend.position='none')+
      geom_segment(aes(x=1, y=2.3, xend=2, yend=2.3))+
      annotate("text", x=1.5, y=2.33, label="**",size=5)+
      annotate("text", x=0.90, y=label, label=
                 sum(n$Var1 %in% pair.filter$pseu))+
      annotate("text", x=1.16, y=label, label=
                 paste0("(",round(sum(n$Var1 %in% pair.filter$pseu)/nrow(n)*100,2),"%)"))+
      annotate("text", x=1.93, y=label, label=
                 length(unique(pair.filter$pseu))-sum(n$Var1 %in% pair.filter$pseu))+
      annotate("text", x=2.16, y=label, label=
                 paste0("(",round(r,2),"%)"))
  }
 
  ggsave(p,filename = file.path("~/Pseudo/Result",S,"Picture",paste0(S,".coexpressedpair.dynamic.pdf")),
         device = "pdf",width = 5,height = 5)
  
}

  
# Enrichment of concomitant dynamics
table(pair.filter$ddpseu,pair.filter$ddcoding>0) %>% fisher.test() #odds ratio = 5.6, p-value <2e-16; odds ratio =173
# while weak dynamic tissue number
cor.test(pair.filter$freq,pair.filter$ddcoding) #cor=0.05, p-value <2e-16; cor=0.07, p-value <2e-16


#####################
# Venn plot
library(VennDiagram)
setwd(file.path("/home/qians/Pseudo/Result",S,"Picture"))
pdf(file.path("~/Pseudo/Result",S,"Picture",paste0(S,".YesorNoDynamic.cod.venn.pdf")))
grid.draw(venn.diagram(list("Dynamic" = unique(pair.filter[pair.filter$ddpseu==1,2]),
                            "Non-dynamic" = unique(pair.filter[pair.filter$ddpseu==0,2])),
                       NULL,
                       imagetype ="svg",
                       fill = c("#d95f0d","#fc9272"),
                       cat.pos=c(0,3),
                       cat.cex=c(2,2),
                       cex = 2,
                       width = 4,height = 4))
dev.off()

#####################
# GO analysis
library(org.Hs.eg.db)
library(clusterProfiler)
## coexpressed pair between non-dynamic pseu and coding
if (S == "Human") {
  nondynamic.cod = setdiff(unique(pair.filter[pair.filter$ddpseu==0,2]),
                           unique(pair.filter[pair.filter$ddpseu==1,2]))
  nondynamic.cod.go = enrichGO(nondynamic.cod,
                               OrgDb = org.Hs.eg.db,
                               keyType = "ENSEMBL",
                               ont = "ALL", 
                               qvalueCutoff = 0.05)
  nondynamic.cod.go.df = as.data.frame(nondynamic.cod.go)
  write.table(nondynamic.cod.go.df,file = 
                file.path("/home/qians/Pseudo/Result",S,"Savedata/DDG/GO",paste0(S,".nondynamic.cod.go.df.csv")),
              row.names = FALSE,sep = ",",quote = FALSE)
  
 #nondynamic.cod.go.df = read.csv(file.path("/home/qians/Pseudo/Result",S,"Savedata/DDG/GO",
 #                                          paste0(S,".nondynamic.cod.go.df.csv")),
 #                                header = T,sep = ",", stringsAsFactors = F)# %>% dplyr::filter(.,ONTOLOGY == "BP")
  nondynamic.cod.go.df = nondynamic.cod.go.df[order(-log10(as.numeric(nondynamic.cod.go.df$qvalue))),]
  nondynamic.cod.go.df$Description %<>% factor(.,levels = .)
  ggplot(nondynamic.cod.go.df,aes(Description,-log10(qvalue)))+theme_classic()+
    geom_bar(stat = "identity",position = "dodge",fill="#40a49a")+ylab("-log10(Q VALUE)")+
    theme(axis.text.y = element_text(size=12),axis.title.y = element_text(size=14),
          axis.text.x = element_text(size=12),axis.title.x = element_text(size=14))+
    scale_y_continuous(expand = c(0, 0),limits = c(0,9))+coord_flip()+
    geom_hline(yintercept = -log10(0.05),color="white",linetype = "dashed")+
    geom_hline(yintercept = -log10(0.01),color="white",linetype = "dashed")
  ggsave(filename=file.path("/home/qians/Pseudo/Result",S,"Picture",paste0(S,".nondynamic.cod.go.df.pdf")),
         device = "pdf",width = 4.5,height = 4)
}

## coexpressed pair between dynamic pseu and coding
if (S == "Human") {
  dynamic.cod = setdiff(unique(pair.filter[pair.filter$ddpseu==1,2]),
                        unique(pair.filter[pair.filter$ddpseu==0,2]))
  dynamic.cod.go = enrichGO(dynamic.cod,
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENSEMBL",
                            ont = "ALL", 
                            qvalueCutoff = 0.05)
  dynamic.cod.go.df = as.data.frame(dynamic.cod.go)
  write.table(dynamic.cod.go.df,file = 
                file.path("/home/qians/Pseudo/Result",S,"Savedata/DDG/GO",paste0(S,".dynamic.cod.go.df.csv")),
              row.names = FALSE,sep = ",",quote = FALSE)
  
  #dynamic.cod.go.df = fread(file.path("/home/qians/Pseudo/Result",S,"Savedata/DDG/GO",paste0(S,".dynamic.cod.go.df.csv")),
  #                             header = T,sep = ",", stringsAsFactors = F)# %>% dplyr::filter(.,ONTOLOGY == "BP")
  dynamic.cod.go.df$qvalue = as.numeric(dynamic.cod.go.df$qvalue)
  dynamic.cod.go.df = dynamic.cod.go.df[order(log10(dynamic.cod.go.df$qvalue)),]
  dynamic.cod.go.df$Description %<>% factor(.,levels = rev(.))
  ggplot(dynamic.cod.go.df,aes(Description,-log10(qvalue)))+theme_classic()+
    geom_bar(stat = "identity",position = "dodge",fill="#40a49a")+ylab("-log10(Q VALUE)")+
    theme(axis.text.y = element_text(size=8),axis.title.y = element_text(size=14),
          axis.text.x = element_text(size=12),axis.title.x = element_text(size=14))+
    scale_y_continuous(expand = c(0, 0))+coord_flip()+
    geom_hline(yintercept = -log10(0.05),color="white",linetype = "dashed")+
    geom_hline(yintercept = -log10(0.01),color="white",linetype = "dashed")+
    scale_x_discrete(labels=function(x) str_wrap(x, width=60))
  ggsave(filename=file.path("/home/qians/Pseudo/Result",S,"Picture",paste0(S,".dynamic.cod.go.df.pdf")),
         device = "pdf",width = 7,height = 11)
}

###merge two above
M = rbind(nondynamic.cod.go.df[1:10,],dynamic.cod.go.df[1:10,])
M = rbind(nondynamic.cod.go.df[order(nondynamic.cod.go.df$qvalue),][1:10,],
          dynamic.cod.go.df[order(dynamic.cod.go.df$qvalue),][1:10,])
M = M[,c(3,4,8,10)]
M$type = c(rep("Nondynamic",min(nrow(nondynamic.cod.go.df),10)),rep("Dynamic",10))

M$All = lapply(M$GeneRatio,function(x)strsplit(x,split = "/",fixed = TRUE)[[1]][2]) %>% unlist() %>% as.numeric()
M$Ratio = M$Count / M$All
M$Description %<>% capitalize()
library(wesanderson)
ggplot(M,aes(type,Description,size = Ratio,color=-log10(qvalue)))+geom_point()+
  #scale_color_manual(values=wes_palette(n=3,type =  "continuous", name="GrandBudapest1"))+
  scale_color_gradient(low="blue", high="red") + theme_classic()+
  theme(axis.text.x = element_text(size = 14, hjust = 1, vjust = 1, angle = 45),
        axis.text.y = element_text(size = 14),
        axis.title = element_blank())+
  #theme(legend.position = c(-2,0.7),legend.background = element_blank())+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))
ggsave(filename=file.path("/home/qians/Pseudo/Result",S,"Picture",paste0(S,".Merge.cod.go.df.pdf")),
       device = "pdf",width = 7.5,height = 5.5)

##simplifyenrichment

a = read.csv(file.path("/home/qians/Pseudo/Result",S,"Savedata/DDG/GO",paste0(S,".dynamic.cod.go.df.xls")),sep = "\t")
    go_id = a$ID
    mat = GO_similarity(go_id,ont = "BP")
    pdf(file = paste0("~/test.pdf"),width = 10,height = 7)
    df = simplifyGO(mat,exclude_words=c("process"))
    dev.off()



## coexpressed pair between common coding
if (TRUE) {
  common.cod = intersect(unique(pair.filter[pair.filter$ddpseu==1,2]),
                        unique(pair.filter[pair.filter$ddpseu==0,2]))
  common.cod.go = enrichGO(common.cod,
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENSEMBL",
                            ont = "ALL", 
                            qvalueCutoff = 0.05)
  common.cod.go.df = as.data.frame(common.cod.go)
  write.table(common.cod.go.df,file = 
                file.path("/home/qians/Pseudo/Result",S,"Savedata/DDG/GO",paste0(S,".common.cod.go.df.csv")),
              row.names = FALSE,sep = ",",quote = FALSE)
  
  common.cod.go.df = read.csv(file.path("/home/qians/Pseudo/Result",S,"Savedata/DDG/GO",
                                         paste0(S,".common.cod.go.df.csv")),
                               header = T,sep = ",", stringsAsFactors = F) %>% 
    dplyr::filter(.,ONTOLOGY == "BP")
  common.cod.go.df$qvalue = as.numeric(common.cod.go.df$qvalue)
  common.cod.go.df = common.cod.go.df[order(-log10(common.cod.go.df$qvalue)),]
  common.cod.go.df$Description %<>% factor(.,levels = .)
  ggplot(common.cod.go.df,aes(Description,-log10(qvalue)))+
    geom_bar(stat = "identity",position = "dodge")+ylab("-log10(Q VALUE)")+
    theme(axis.text.y = element_text(size=18),
          axis.text.x = element_text(size=18))+theme_classic()+
    scale_y_continuous(expand = c(0, 0))+coord_flip()+
    geom_hline(yintercept = -log10(0.05),color="white",linetype = "dashed")+
    geom_hline(yintercept = -log10(0.01),color="white",linetype = "dashed")+
    scale_x_discrete(labels=function(x) str_wrap(x, width=60))
  ggsave(filename=file.path("/home/qians/Pseudo/Result",S,"Picture",paste0(S,".common.cod.go.df.pdf")),
         device = "pdf",width = 7,height = 10)
}


  
#####################
#DOSE analysis
#https://songqi.netlify.app/2020/04/18/go_kegg_do/





#####################
#cancer

rm(pair)
cancer = read.csv("/home/qians/Pseudo/Data/Ref/Human/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv",
                  header = T,sep = "\t",stringsAsFactors = F)
library(org.Hs.eg.db)
sym = toTable(org.Hs.egSYMBOL)
ens = toTable(org.Hs.egENSEMBL)
sym$ens = ens[match(sym$gene_id,ens$gene_id),2]
cancer$ens = sym[match(cancer$SYMBOL,sym$symbol),3]
cancer[cancer$SYMBOL=="FAM46C","ens"]="ENSG00000183508"

# Cancer genes have less co-expressed pair
freq = table(pair.filter$coding) %>% as.data.frame()
freq$cancer = cancer[match(freq$Var1,cancer$ens),1]
freq[is.na(freq)] = "No"
freq[freq$cancer!="No","cancer"]= "Yes"
ggplot(freq,aes(cancer,Freq))+geom_boxplot()
# dynamic pseudogenes show no bias to co-express with cancer gene
pair.filter$cancer = cancer[match(pair.filter$coding,cancer$ens),1]
pair.filter[is.na(pair.filter)] = 0
pair.filter[pair.filter$cancer!=0,"cancer"]=1
table(pair.filter$cancer,pair.filter$ddpseu)

pair.filter.cancer = na.omit(pair.filter)

######################
#constraint 
#We characterized DDGs using three different metrics of functional constraint: (1) the residual variation intolerance score (RVIS); (2) the probability of being intolerant to loss-of-function mutations (pLI score); and (3) the selection against heterozygous loss of gene function (shet). All metrics were applied to data from the Exome Aggregation Consortium (ExAC)9. We obtained the pLI and RVIS scores from ref. 7 and shet values from ref. 10. We also used the copy number variation (CNV) intolerance score as applied to the ExAC data from ref. 11. The lists of transcription factors were from the animalTFDB (v.2.0)52.
#Gene expression across mammalian organ development
#High-throughput discovery of novel developmental phenotypes
#Genic Intolerance to Functional Variation and the Interpretation of Personal Genomes


