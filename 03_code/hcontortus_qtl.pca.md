



## PCA of nuclear and mitochondrial variants

```bash
cd /lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL/05_ANALYSIS/PCA

# mtDNA
ln -s ../../04_VARIANTS/FILTERED/HCON_QTL.cohort.2023-12-12.n278.mito_variants.final.recode.vcf mito.vcf

```

vcftools --vcf mito.vcf --missing-indv

cat out.imiss | awk '{if($5<0.01) print $1}' > keep.list

vcftools --gzvcf REFsplit.7.list.cohort.tmp.vcf.gz --keep keep.list --max-missing 1 --remove-indels --recode --out hcontortus_chr_mtDNA_arrow_pilon.missindv0.1.maxmiss1.n278

#After filtering, kept 278 out of 590 Individuals
#After filtering, kept 1033 out of a possible 3403 Sites


grep "XQTL\|MHCO3\|GB_ISE" keep.list > keep.XQTL.list

vcftools --gzvcf REFsplit.7.list.cohort.tmp.vcf.gz  --keep keep.XQTL.list --max-missing 1 --remove-indels --recode --out hcontortus_chr_mtDNA_arrow_pilon.missindv0.1.maxmiss1.n256.xqtl

After filtering, kept 256 out of 590 Individuals
After filtering, kept 1250 out of a possible 3403 Sites

```R
library(tidyverse)
library(SNPRelate)
library(ggrepel)


#vcf.fn <- "hcontortus_chr_mtDNA_arrow_pilon.raw.vcf.gz"
#vcf.fn <-  "hcontortus_chr_mtDNA_arrow_pilon.missindv0.1.maxmiss1.n276.recode.vcf"
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


ggplot(data, aes(EV1, EV2, colour=population, label=sample.id)) + 
     geom_point() +
     theme_bw() +
     geom_text_repel(data = subset(data, population == "Parent_susceptible"), max.overlaps = Inf) +
     labs(title="mitochondrial_variants",
          x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
          y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"))

```

## nuclear variants
```R
library(tidyverse)
library(SNPRelate)

vcf.fn <- "hcontortus_chr1_Celeg_TT_arrow_pilon.raw.vcf.gz"

snpgdsVCF2GDS(vcf.fn, "nuclear.gds", method="biallelic.only")

genofile <- snpgdsOpen("nuclear.gds")

pca <- snpgdsPCA(genofile, num.thread=2, autosome.only=F)

pca <- snpgdsPCA(genofile, num.thread=2, autosome.only=F)

data <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(data)


population <- read.table("sample-id_populations.txt", header=T, sep="\t") 

data <- inner_join(data, population, by="sample.id")


ggplot(data, aes(EV1, EV2, colour=population)) + 
     geom_point() +
     labs(title="mitochondrial_variants",
          x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
          y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"))
```