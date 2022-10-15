rm(list = ls());gc();rm(list = ls())
Num = "015.1."

library(ImmuLncRNA)
mRNA_exp_data <- system.file("extdata","mRNA_test.txt",package = "ImmuLncRNA")
mRNA_exp_data

library(estimate)
dir_test <- "~/R/x86_64-pc-linux-gnu-library/4.0/ImmuLncRNA/extdata"
file_test <-"mRNA_test.txt"
res_test <- ImmuLncRNA::turpur.est(dir_test,file_test)
res_test[1:3] # showing this 3 samples' tumor purity.

lncRNA_exp <- as.matrix(lncRNA_exp)
mRNA_exp <- as.matrix(mRNA_exp)
turpur_ori <- as.numeric(turpur[,1])
names(turpur_ori) <- rownames(turpur)
lncRNA_exp[1:3,1:2]
#>                 TCGA-OR-A5LE-01A-11R-A29S-07 TCGA-OR-A5JM-01A-11R-A29S-07
#> ENSG00000082929                    14.396825                     9.696367
#> ENSG00000093100                     8.332295                    10.386720
#> ENSG00000099869                    19.338013                    17.315262
mRNA_exp[1:3,1:2]
#>                 TCGA-OR-A5LE-01A-11R-A29S-07 TCGA-OR-A5JM-01A-11R-A29S-07
#> ENSG00000000003                    18.259242                     17.80044
#> ENSG00000000005                     8.758933                     10.22840
#> ENSG00000000419                    19.585016                     19.85835
turpur_ori[1:2]
#> TCGA-OR-A5LE-01A-11R-A29S-07 TCGA-OR-A5JM-01A-11R-A29S-07 
#>                    0.9945247                    0.9883493

test_res <- par.cor(mRNA_exp,lncRNA_exp,turpur_ori,adjusted=T)
class(test_res)
#> [1] "list"
test_res$pcor.value[1:3,1:3]
#>                 ENSG00000000003 ENSG00000000005 ENSG00000000419
#> ENSG00000082929      0.11242325     -0.11147442       0.1584019
#> ENSG00000093100      0.02499571     -0.20038392      -0.2302502
#> ENSG00000099869      0.01908235      0.08885238      -0.1390293
test_res$p.value[1:3,1:3]
#>                 ENSG00000000003 ENSG00000000005 ENSG00000000419
#> ENSG00000082929       0.3239682      0.32812060      0.16194095
#> ENSG00000093100       0.8274481      0.07457463      0.03914147
#> ENSG00000099869       0.8678529      0.43676185      0.22098407

lncRNA_exp <- as.matrix(lncRNA_exp)
mRNA_exp <- as.matrix(mRNA_exp)
turpur_ori <- as.numeric(turpur[,1])
names(turpur_ori) <- rownames(turpur)
# str(pathways) : check for pathway list
k=0.995
library(ImmuLncRNA)
test_res <- immu.LncRNA(mRNA_exp,lncRNA_exp,adjusted=T,turpur_ori,pathways,k)
#> Loading required package: Rcpp
test_res$sig_pairs[1:2,] # showing the immu-lncRNA pairs
#>      lncRNA            pathway     
#> [1,] "ENSG00000082929" "Chemokines"
#> [2,] "ENSG00000082929" "Cytokines"
test_res$fgseaRes_all[1:2,] # showing the GSEA results
#>             lncRNA                             pathway      pval padj
#> 1: ENSG00000082929 Antigen_Processing_and_Presentation 0.3867735    1
#> 2: ENSG00000082929                      Antimicrobials 0.6032064    1
#>            ES        NES nMoreExtreme size     leadingEdge   sigValue
#> 1: -0.8421053 -1.1054384          192    1 ENSG00000001167 -0.2264529
#> 2: -0.7368421 -0.9672586          300    1 ENSG00000000938  0.2064128


