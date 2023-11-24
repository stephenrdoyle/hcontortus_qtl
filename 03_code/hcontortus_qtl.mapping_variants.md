# hcontortus_qtl

The primary focus of this project on QTL mapping ivermectin response in Haemonchus contortus.

The origin of the material analysed is the MHco3/MHco18 genetic cross performed by the BUG Consortium. The analysis of these data have been described previously in different ways, including
- genetic map (Doyle et al. 2018 GBE)
- XQTL genomics (Doyle, Laing et al. 2022 Cell Reports)
- XQTL transcriptomics (Laing et al. 2022 PLoS Pathogens)

The main QTL experiment here is a dose response larval development experiment, in which two groups of larvae were collected:
- L3 that developed normally at high doses of ivermectin
- L1/L2 that were not developing well on low doses of ivermectin
We originally performed a poolseq experiment, where pools of larvae (n=200 per pool) were sequenced and compared, revealing a consistent chromosome 5 signature of selection as the XQTL experiments. The rationale for this new QTL experiment was that genotypes from individual worms, rather than allele frequencies from pooled worms, may provide additional support for variants in the region under selection.

The added benefit of sequencing individual larvae is that there a lot more than can be done with the data, including (and this forms my to-do list)
- determine LD and explore recombination
- perform phasing
- sex-specific analyses, which can be determined from X chromosome coverage and heterozygosity
- check for variation in ploidy
- genetic relatedness
- more accurate CNV analyses
- perform GWAS, taking into account genetic relatedness.

Further, this is going to be an excellent resource for some benchmarking, including
- SNP calling
- genome graphs


We also prepared sequencing libraries from larvae obtained from US farms that have been phenotypically tested using drenchrite assays by Ray Kaplan. These are the same farms for which we have sequenced pools and analysed in the XQTL paper. Library prep revealed these didn't work that well - the libraries were quite weak especially compared with the QTL samples - and so we have taken a different pooling approached to try and maximise samples. There were more samples prepared than sequenced, however, becasue of the library strength we have put these on hold for now. For those samples that were sequecned, these probably wont work that well.


## Sequencing data  
### Sequencing data - QTL
- susceptible (267 samples)
     - 36342_3 - 94 samples (SUS)
     - 36808_1 - 79 samples (SUS) - QC PENDING
     - 35990_1 - 94 samples (SUS)
- resistance (235 samples)
     - 36342_1 - 75 samples (RES)
     - 36342_2 - 82 samples (RES)
     - 36342_4 - 78 samples (RES)


### Sequencing data - Farm
- 35887_2 - 96 samples (FARM)




```bash
# get reference
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/QTL/01_REFERENCE

cp ../../GENOME/REF/HAEM_V4_final.chr.fa .

# make a index and a dict file
samtools faidx HAEM_V4_final.chr.fa
samtools dict HAEM_V4_final.chr.fa > HAEM_V4_final.chr.dict
```
[↥ **Back to top**](#top)



```bash
# get the raw data
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/QTL/02_RAW

for i in 36342_3 36808_1 35990_1 36342_1 36342_2 36342_4 35887_2; do
     pf data --type lane --id ${i} --symlink ./ --rename --filetype fastq;
     done

# extract lane/barcode and sample IDs to create some metadata
>sample_lanes.list
for i in 36342_3 36808_1 35990_1 36342_1 36342_2 36342_4 35887_2; do
     pf supplementary --type lane --id ${i} | grep -v "Sample" | sed 's/\#/_/g' | awk '{print $3,$6}' OFS="\t" >> sample_lanes.list;
     done

```
[↥ **Back to top**](#top)



```bash
# mapping
cd /nfs/users/nfs_s/sd21/lustre118_link/hc/QTL/03_MAPPING

ln -s ../02_RAW/sample_lanes.list

screen
while read LANE NAME; do
     ~sd21/bash_scripts/run_bwamem_splitter \
     ${NAME} \
     /nfs/users/nfs_s/sd21/lustre118_link/hc/XQTL/01_REFERENCE/HAEM_V4_final.chr.fa \
     /nfs/users/nfs_s/sd21/lustre118_link/hc/QTL/02_RAW/${LANE}_1.fastq.gz \
     /nfs/users/nfs_s/sd21/lustre118_link/hc/QTL/02_RAW/${LANE}_2.fastq.gz;
     done < sample_lanes.list &


```





## merge UGA and ISE split mapping into single bams

```bash
ls -1d US_UGA* | cut -c-10 | sort | uniq > samples.list && ls -1d GB_ISE* | cut -c-10 | sort | uniq >> samples.list

while read -r NAME; do 
mkdir ${NAME}_merged
ls -1 $PWD/${NAME}.[1-4]/${NAME}.[1-4].bam > ${NAME}_merged/bam.list
samtools merge -o ${NAME}_merged/${NAME}.merged.bam -b ${NAME}_merged/bam.list;
samtools sort ${NAME}_merged/${NAME}.merged.bam -o ${NAME}_merged/${NAME}.sorted.bam; 
samtools index ${NAME}_merged/${NAME}.sorted.bam;
rm ${NAME}_merged/${NAME}.merged.bam;
samtools flagstat ${NAME}_merged/${NAME}.sorted.bam > ${NAME}_merged/${NAME}.sorted.flagstat;
done < samples.list

# note: flagstat only useful for total reads, as only good quality reads were kept in the intial mapping
```




## extract Kraken data

```bash
cd /lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL/05_ANALYSIS/KRAKEN 

cp ../../lanes_samples.list .

# pull kraken qc reports
cat lanes_samples.list | cut -c-7 | sort | uniq | while read -r LANE; do
     pf qc --type lane --id ${LANE} --symlink ./ --rename; 
done

# fix the lanes and samples file
sed -i 's/_/#/2' lanes_samples.list

# copy the symlinked file and change the name of the file from lane to sample name
cat lanes_samples.list | while read -r LANE NAME; do
     cp ${LANE}_kraken.report ${NAME}_kraken.report;
done

# fix the header by removing the first two lines - only a fix for Sanger kraken reports
for i in *_kraken.report; do 
     sed -i '1,2d' ${i}; 
     done

# remove the symlinked files
rm 3*

# make a multiqc report of the kraken data
multiqc .
```






## PCA of nuclear and mitochondrial variants

```bash
cd /lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL/05_ANALYSIS/PCA

# mtDNA
ln -s ../../04_VARIANTS/gatk_hc_test/GATK_HC_MERGED/hcontortus_chr_mtDNA_arrow_pilon.raw.vcf.gz

```



```R
library(tidyverse)
library(SNPRelate)


vcf.fn <- "hcontortus_chr_mtDNA_arrow_pilon.raw.vcf.gz"

snpgdsVCF2GDS(vcf.fn, "mtDNA.gds", method="biallelic.only")

snpgdsSummary("mtDNA.gds")

genofile <- snpgdsOpen("mtDNA.gds")

pca <- snpgdsPCA(genofile, num.thread=2, autosome.only=F)



pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

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