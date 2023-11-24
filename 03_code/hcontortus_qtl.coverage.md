# H.contortus QTL: Coverage analyses


```bash
# make links to all the bams
cd /lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL

OUTPUT_DIR=/lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL/05_ANALYSIS/COVERAGE
find ~+ -type f -name '*.bam*' -exec ln -vs "{}" $OUTPUT_DIR/ ';'
rm BLANK*

# switch directories
cd /lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL/05_ANALYSIS/COVERAGE

# get Hcon annotation from WBPS
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3.gz


#generate the bed file for exons
zcat haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3.gz | grep "exon" | awk -F'[:\t]' '{print $1, $4, $5}' OFS="\t" > exon.bed

for i in *.bam; do
    samtools bedcov -Q 20 exon.bed ${i} | awk -F'\t' '{printf $1"\t"$2"\t"$3"\t"$4"\t"$4/($3-$2)"\n"}' OFS="\t" > ${i%.bam}.exon.cov
done

# extract mtDNA and nuclear (median & stddev) data. For exon data.
for i in *.exon.cov; do
     name=${i%.exon.cov};
     autosome=$(grep "hcontortus_chr1_Celeg_TT_arrow_pilon\|hcontortus_chr2_Celeg_TT_arrow_pilon\|hcontortus_chr3_Celeg_TT_arrow_pilon\hcontortus_chr4_Celeg_TT_arrow_pilon\|hcontortus_chr5_Celeg_TT_arrow_pilon)
     ]" ${i} | datamash median 5 sstdev 5 );
     xchromosome=$(grep "hcontortus_chrX_Celeg_TT_arrow_pilon" ${i} | datamash median 5 sstdev 5);
     echo -e "${name}\t${autosome}\t${xchromosome}";
done > autosome_to_Xchromsome_cov_exon.stat

```



```R

#plotting the data median and reads mapped all samples
library(tidyverse)
cov_data <- read.table("autosome_to_Xchromsome_cov_200kb_med.stat")
meta_data <- read.table ("all_samples_metadata2.txt")
data <- left_join(cov_data,meta_data,by="V1")
colnames(data) <- c("sample_name", "autosome_cov_median", "autosome_cov_sd", "x_cov_median", "x_cov_sd", "reads_mapped")
ggplot(data, aes(sample_name, x_cov_median/autosome_cov_median)) +
      geom_point(aes(colour = log10(as.numeric(reads_mapped)))) +
      scale_color_viridis_c()+
      labs(title="X-to-autosomal coverage ratio all samples 100kb median", x="Sample name", y="X-to-autosomal coverage ratio") +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      theme_bw(base_size = 10) +
      ylim(0.1,2) +
      coord_flip()
ggsave("plot_x-to-autosome_ratio_sexdet_all.png")
ggsave("plot_x-to-autosome_ratio_sexdet_all.pdf",height=10, width=7, useDingbats=FALSE)




```