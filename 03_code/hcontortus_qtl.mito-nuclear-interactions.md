# Hcontortus_qtl: mito-nuclear-interactions

### Stephen Doyle

## background
- PCA analysis of the QTL cohort has reveal at least 3 groups of samples based on mitochondrial variants, suggestign segregating haplotypes from the ISE parental population used to intitiate the cross
- however, the nuclear variants are admixed - this is expected
- I want to now explore if there are any interesting differences in the nuclear genome when comparing samples based on the mtDNA haplogroup
    - expect no major differences
    - however, if there are, it may suggest there is an interaction between mtDNA and nuclear genomes



## working directory
```bash
cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/QTL/05_ANALYSIS/PARENT_BIAS/MITO-NUCLEAR-INTERACTION

```


## defining the sample sets based on mtDNA groups
- first need to recreate the PCA, and then create lists of sample IDs based on 
```bash
# recreate the mtDNA PCA

ln -s ../PCA/hcontortus_chr_mtDNA_arrow_pilon.missindv0.1.maxmiss1.n256.xqtl.recode.vcf


```

```R
R
#R version 4.0.3

library(tidyverse)
library(SNPRelate)
library(ggrepel)


vcf.fn <- "hcontortus_chr_mtDNA_arrow_pilon.missindv0.1.maxmiss1.n256.xqtl.recode.vcf" 

snpgdsVCF2GDS(vcf.fn, "mtDNA_xqtl.gds", method="biallelic.only")

snpgdsClose(genofile)
genofile <- snpgdsOpen("mtDNA_xqtl.gds")

pca <- snpgdsPCA(genofile, num.thread=2, autosome.only=F)

data <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)


population <- read.table("sample-id_populations.txt", header=T, sep="\t") 

data <- inner_join(data, population, by="sample.id")


plot <- ggplot(data, aes(EV1, EV2, label=sample.id)) + 
     geom_point() +
     theme_bw() +
     labs(title="groups of mitochondrial_variants",
          x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
          y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"))

plot

# add some rectangles to the plot to define groups
plot + 
# top left cluster
    geom_rect(aes(xmin=-0.125,xmax=-0.05, ymin=0 , ymax=0.05), fill=NA, col="grey") + 
    geom_text(aes(x=-0.125,y=0.05, label="group 1"), hjust = 0) +
# top right cluster
    geom_rect(aes(xmin=0.03, xmax=0.085, ymin=0.025, ymax=0.075), fill=NA, col="grey") + 
    geom_text(aes(x=0.03,y=0.075, label="group 2"), hjust = 0) +
# middle cluster
    geom_rect(aes(xmin=-0.03, xmax=0.025, ymin=-0.04, ymax=0.01), fill=NA, col="grey") + 
    geom_text(aes(x=-0.03,y=0.01, label="group 3"), hjust = 0) +
# bottom cluster
    geom_rect(aes(xmin=-0.01, xmax=0.05, ymin=-0.17, ymax=-0.06), fill=NA, col="grey") + 
    geom_text(aes(x=-0.01,y=-0.06, label="group 4"), hjust = 0)

ggsave("mtDNA.pca.mito-groups.png")
```
![](../04_analysis/mtDNA.pca.mito-groups.png)
