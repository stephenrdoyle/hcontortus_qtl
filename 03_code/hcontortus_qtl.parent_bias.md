# H.contortus QTL: Parental bias 

### Author: Stephen Doyle

- the XQTL progeny sequenced is the F5 generation of the cross between MHco3 and MHco18
- therefore, we'd expect that by this time point, after multiple rounds of passage, mating, and recombination, that the genomes would be sufficiently admixed
- however, there might be parts of the genome, that, for whatever reason might look more like one or the other parents, if there is some degree of "hidden" selection going on
    - this this is the case, it would be interesting to know what genes might be being selected for
    - note this should be independent of drug resistance from MHco18, given these progeny have not seen drug until the F5 generation
        - the XQTL_resistant samples will have a biased distribution of MHco18 alleles associated with IVM resistance
        - however, the XQTL_susceptible samples should not


### Aim
- to determine if there is parental allele bias
- to determine if there is mito-nuclear discordance based on mitochondrial haplotype


### Approach
- identify variants that are fixed between parental strains, and then determine if they are segregating as expected in the cross
    - as expected might mean intermediate allele frequencies
    - it might also mean stable under hardy weinberg equilibrium
- for the haplotype analysis, scan genome-wide between groups based on the mtDNA grouping.



## Identify variants with an allelic bias between the parental MHco3 and MHco18 strains
### Calculate allele frequencies of the parental strains
```bash 
# working dir 
cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/QTL/05_ANALYSIS/PARENT_BIAS

# get BAMs of the parental strains - these are from the XQTL paper (Doyle et al 2022 Cell Reports) and derived from pools of 200 larvae
ln -s ../../QTL_MAPPING/MHCO3_P0_L3_n200_01/MHCO3_P0_L3_n200_01.bam
ln -s ../../QTL_MAPPING/MHCO3_P0_L3_n200_01/MHCO3_P0_L3_n200_01.bam.bai
ln -s ../../QTL_MAPPING/MHCO18_P0_L3_n200_IVM_01/MHCO18_P0_L3_n200_IVM_01.bam
ln -s ../../QTL_MAPPING/MHCO18_P0_L3_n200_IVM_01/MHCO18_P0_L3_n200_IVM_01.bam.bai



# calculate variant frequency and coverage of the two parental strains
#-- using grenedalf, whcih is great new tool for poolseq data, and replacement for popoolation2

module load grenedalf/0.2.0

bsub.py --queue long 10 grenedalf_freq \
"grenedalf frequency \
    --write-sample-coverage \
    --write-sample-alt-freq \
    --write-total-frequency \
    --sam-min-map-qual 30 \
    --sam-min-base-qual 30 \
    --file-prefix MHCO3_v_MHCO18_pools \
    --separator-char tab \
    --sam-path MHCO3_P0_L3_n200_01.bam \
    --sam-path MHCO18_P0_L3_n200_IVM_01.bam"


# once completed, calculate median coverage of each group
#--- want this data, as coverage will be variable, and allele frequencies will be somewhat biased by low coverage
# total
cat MHCO3_v_MHCO18_poolsfrequency.csv | datamash --header-out --header-in median 5,7

#median(MHCO18_P0_L3_n200_IVM_01:1.COV)	median(MHCO3_P0_L3_n200_01:1.COV)
#58	69

# by chromosome
cat MHCO3_v_MHCO18_poolsfrequency.csv | datamash --group 1 --header-out --header-in median 5,7

#GroupBy(CHROM)	median(MHCO18_P0_L3_n200_IVM_01:1.COV)	median(MHCO3_P0_L3_n200_01:1.COV)
#hcontortus_chr1_Celeg_TT_arrow_pilon	55	66
#hcontortus_chr2_Celeg_TT_arrow_pilon	63	77
#hcontortus_chr3_Celeg_TT_arrow_pilon	56	66
#hcontortus_chr4_Celeg_TT_arrow_pilon	54	68
#hcontortus_chr5_Celeg_TT_arrow_pilon	54	68
#hcontortus_chrX_Celeg_TT_arrow_pilon	68	70
#hcontortus_chr_mtDNA_arrow_pilon	1609	1587
```

### Identify variants that polarise the two strains
- need to decide
    - a sensible allele frequency difference - decided on a difference of 0.8
    - a minimum coverage - decided on 1/3 of the median read depth

```bash
# checking variants with biased allele frequency
#-- freq diff = 0.8
#-- min coverage = 0.33 * median
awk '{print $1,$2,$8,$7,$6,$5,$6-$8}' OFS="\t" MHCO3_v_MHCO18_poolsfrequency.csv |\
    awk '{if($7>0.8 && $4>23 && $6>19) print}' OFS="\t" |\
    cut -f1 | sort | uniq -c

# #CHR POS MHCO3_freq MHCO3_cov MHCO18_freq MHCO18_cov diff(MHCO18_freq-MHCO3_freq)

# 0.8 cutoff - 35711 total
   7024 hcontortus_chr1_Celeg_TT_arrow_pilon
  12695 hcontortus_chr2_Celeg_TT_arrow_pilon
   2247 hcontortus_chr3_Celeg_TT_arrow_pilon
   4808 hcontortus_chr4_Celeg_TT_arrow_pilon
   8244 hcontortus_chr5_Celeg_TT_arrow_pilon
    693 hcontortus_chrX_Celeg_TT_arrow_pilon 

# also tested a 0.9 cutoff, to check variant numbers
awk '{print $1,$2,$8,$7,$6,$5,$6-$8}' OFS="\t" MHCO3_v_MHCO18_poolsfrequency.csv |\
    awk '{if($7>0.9 && $4>23 && $6>19) print}' OFS="\t" |\
    cut -f1 | sort | uniq -c

# 0.9 cutoff - 10706 total
   2166 hcontortus_chr1_Celeg_TT_arrow_pilon
   3721 hcontortus_chr2_Celeg_TT_arrow_pilon
    568 hcontortus_chr3_Celeg_TT_arrow_pilon
   1282 hcontortus_chr4_Celeg_TT_arrow_pilon
   2811 hcontortus_chr5_Celeg_TT_arrow_pilon
    158 hcontortus_chrX_Celeg_TT_arrow_pilon



# checking what the distribution of SNPs look like - extract the postions for plotting
awk '{print $1,$2,$8,$7,$6,$5,$6-$8}' OFS="\t" MHCO3_v_MHCO18_poolsfrequency.csv |\
    awk '{if($7>0.8 && $4>23 && $6>19) print $1,$2,$7,"UGA"; else if ($7<-0.8 && $4>23 && $6>19) print $1,$2,$7,"ISE"}' OFS="\t" > MHCO3_v_MHCO18_poolsfrequency.0.8freq.pos

awk '{print $1,$2,$8,$7,$6,$5,$6-$8}' OFS="\t" MHCO3_v_MHCO18_poolsfrequency.csv |\
    awk '{if($7>0.9 && $4>23 && $6>19) print $1,$2,$7,"UGA"; else if ($7<-0.9 && $4>23 && $6>19) print $1,$2,$7,"ISE"}' OFS="\t" > MHCO3_v_MHCO18_poolsfrequency.0.9freq.pos
```


### make a plot
```R
library(tidyverse)


data_0.8 <- read.table("MHCO3_v_MHCO18_poolsfrequency.0.8freq.pos", header=F)

plot_0.8 <- ggplot(data_0.8, aes(V2*1e6,abs(V3),col=V1)) + 
    geom_point(size=0.1) + 
    facet_grid(V1~V4) + 
    ylim(0,1) +
    theme_bw() +
    guides(color = FALSE) +
    labs(title="Variant positions that differ between MHco3 and MHco18 by > 0.8", 
        y="Variant frequency difference (min: MHCO18_freq-MHCO3_freq >0.8)", 
        x="Chromosomal position")
plot_0.8

ggsave("MHCO3_v_MHCO18_poolsfrequency.0.8freq.pos.pdf", height=7, width=7)
ggsave("MHCO3_v_MHCO18_poolsfrequency.0.8freq.pos.png")

data_0.9 <- read.table("MHCO3_v_MHCO18_poolsfrequency.0.9freq.pos", header=F)

plot_0.9 <- ggplot(data_0.9, aes(V2*1e6,V3,col=V1)) + 
    geom_point(size=0.1) + 
    facet_grid(V1~.) + 
    ylim(0,1) +
    theme_bw() +
    guides(color = FALSE) +
    labs(title="Variant positions that differ between MHco3 and MHco18 by > 0.9", 
        y="Variant frequency difference (min: MHCO18_freq-MHCO3_freq >0.9)", 
        x="Chromosomal position")
plot_0.9

ggsave("MHCO3_v_MHCO18_poolsfrequency.0.9freq.pos.pdf", height=7, width=7)
ggsave("MHCO3_v_MHCO18_poolsfrequency.0.9freq.pos.png")
```
![](../04_analysis/MHCO3_v_MHCO18_poolsfrequency.0.8freq.pos.png)
![](../04_analysis/MHCO3_v_MHCO18_poolsfrequency.0.9freq.pos.png)

- from the plots, it is clear that these variants that discriminate MHco3 and MHco18 are not randomly distributed
    - knowing where the QTL for BZ, LEV, and IVM are from the XQTL paper (Doyle 2022 Cell Rep, it is clear from these data there are enriched regions of "fixed" alleles
        - IVM - 37.5 on chr 5, possibly start of chr2
        - BZ - 7 Mb on chr 1 and 14 Mb on Chr2
        - LEV - 30 Mb on chr5
    - however, there are other regions clearly that are different
    - also regions where there are suspicionsly few variants
        - X chromosome a good example - note, might be more impacted by coverage cutoffs

- will stick with the 0.33X coverage cutoffs and diff of 0.8


```bash
# extract positions to keep for filtering the VCF containing all samples

awk '{print $1,$2,$8,$7,$6,$5,$6-$8}' OFS="\t" MHCO3_v_MHCO18_poolsfrequency.csv | awk '{if($7>0.8 && $4>23 && $6>19) print $1,$2}' OFS="\t" > MHCO3_v_MHCO18_diff0.8_biased_variants.keep-positions

awk '{print $1,$2,$8,$7,$6,$5,$6-$8}' OFS="\t" MHCO3_v_MHCO18_poolsfrequency.csv | awk '{if($7>0.8 && $4>23 && $6>19) print $1,$2 ; else if ($7<-0.8 && $4>23 && $6>19) print $1,$2}' OFS="\t" > MHCO3_v_MHCO18_diff0.8_biased_variants.keep-positions

```





########  TESTING  #########

ln -s ../../04_VARIANTS/gatk_hc_test/GATK_HC_MERGED/hcontortus_chr5_Celeg_TT_arrow_pilon.raw.vcf.gz
ln -s ../../04_VARIANTS/gatk_hc_test/GATK_HC_MERGED/hcontortus_chr1_Celeg_TT_arrow_pilon.raw.vcf.gz

XQTL_resistant.samples.keep_list
XQTL_susceptible.samples.keep_list - n = 219

# XQTL_susceptible - freq
vcftools \
    --gzvcf hcontortus_chr5_Celeg_TT_arrow_pilon.raw.vcf.gz \
    --keep XQTL_susceptible.samples.keep_list \
    --positions MHCO3_v_MHCO18_diff0.8_biased_variants.keep-positions \
    --min-alleles 2 \
    --max-alleles 2 \
    --freq \
    --out chr5_XQTL_susceptible_diff0.8_biased_variants

sed -i -e 's/:/\t/g' -e '1d' chr5_XQTL_susceptible_diff0.8_biased_variants.frq


vcftools \
    --gzvcf hcontortus_chr1_Celeg_TT_arrow_pilon.raw.vcf.gz \
    --keep XQTL_susceptible.samples.keep_list \
    --positions MHCO3_v_MHCO18_diff0.8_biased_variants.keep-positions \
    --min-alleles 2 \
    --max-alleles 2 \
    --freq \
    --out chr1_XQTL_susceptible_diff0.8_biased_variants

sed -i -e 's/:/\t/g' -e '1d' chr1_XQTL_susceptible_diff0.8_biased_variants.frq











# XQTL_resistant - freq
vcftools \
    --gzvcf hcontortus_chr5_Celeg_TT_arrow_pilon.raw.vcf.gz \
    --keep XQTL_resistant.samples.keep_list \
    --positions MHCO3_v_MHCO18_diff0.8_biased_variants.keep-positions \
    --min-alleles 2 \
    --max-alleles 2 \
    --freq \
    --out chr5_XQTL_resistant_diff0.8_biased_variants

sed -i -e 's/:/\t/g' -e '1d' chr5_XQTL_resistant_diff0.8_biased_variants.frq


vcftools \
    --gzvcf hcontortus_chr1_Celeg_TT_arrow_pilon.raw.vcf.gz \
    --keep XQTL_resistant.samples.keep_list \
    --positions MHCO3_v_MHCO18_diff0.8_biased_variants.keep-positions \
    --min-alleles 2 \
    --max-alleles 2 \
    --freq \
    --out chr1_XQTL_resistant_diff0.8_biased_variants

sed -i -e 's/:/\t/g' -e '1d' chr1_XQTL_resistant_diff0.8_biased_variants.frq








# XQTL_susceptible - hardy
vcftools \
    --gzvcf hcontortus_chr5_Celeg_TT_arrow_pilon.raw.vcf.gz \
    --keep XQTL_susceptible.samples.keep_list \
    --positions MHCO3_v_MHCO18_diff0.8_biased_variants.keep-positions \
    --min-alleles 2 \
    --max-alleles 2 \
    --hardy \
    --out chr5_XQTL_susceptible_diff0.8_biased_variants

# XQTL_resistant - hardy
vcftools \
    --gzvcf hcontortus_chr5_Celeg_TT_arrow_pilon.raw.vcf.gz \
    --keep XQTL_resistant.samples.keep_list \
    --positions MHCO3_v_MHCO18_diff0.8_biased_variants.keep-positions \
    --min-alleles 2 \
    --max-alleles 2 \
    --012 \
    --out chr5_XQTL_resistant_diff0.8_biased_variants



awk '{$1=""; print $0}' OFS="\t" chr5_XQTL_resistant_diff0.8_biased_variants.012 > chr5_XQTL_resistant_diff0.8_biased_variants.012.2

paste chr5_XQTL_resistant_diff0.8_biased_variants.012.indv chr5_XQTL_resistant_diff0.8_biased_variants.012.2 | awk '{for (i=1; i<=NF; i++) a[i]=a[i](NR!=1?FS:"")$i} END {for (i=1; i in a; i++) print a[i]}' OFS="\t" > 012.indv.geno

rm tmp.012.pos2
printf "CHR\tPOS\n" > tmp.012.pos2 ; cat chr5_XQTL_resistant_diff0.8_biased_variants.012.pos >> tmp.012.pos2

paste tmp.012.pos2 012.indv.geno > genotype_matrix



while read line; do 
coords=$(echo ${line} | cut -f1,2 -d " ")
missing=$(echo ${line} | cut -f3- -d " " | grep -o "\-1" | wc -l)
homozygous_reference=$(echo ${line} | cut -f3- -d " " | grep -o "0" | wc -l)
heterozygous=$(echo ${line} | cut -f3- -d " " | grep -o " 1" | wc -l)
homozygous_variant=$(echo ${line} | cut -f3- -d " " | grep -o "2" | wc -l)

total=$(echo "scale=3; $homozygous_reference + $heterozygous + $homozygous_variant" | bc -l)

#p=$(((2*${homozygous_reference}+${heterozygous})/${total}))
#q=$(((2*${homozygous_variant}+${heterozygous})/${total}))

p=$(echo "scale=3; ((2 * ${homozygous_reference}) + ${heterozygous}) / (2 * ${total})" | bc -l)
q=$(echo "scale=3; ((2 * ${homozygous_variant}) + ${heterozygous}) / (2 * ${total})" | bc -l)

exp_homo_ref=$(echo "scale=3; (${p}^2) * ${total}" | bc -l)
exp_het=$(echo "scale=3; (2*${p}*${q}) * ${total}" | bc -l)
exp_homo_var=$(echo "scale=3; (${q}^2) * ${total}" | bc -l)

chi_hom_ref=$(echo "scale=3; (${homozygous_reference}-${exp_homo_ref})^2 / ${exp_homo_ref}" | bc -l)
chi_het=$(echo "scale=3; (${heterozygous}-${exp_het})^2/${exp_het}" | bc -l)
chi_hom_var=$(echo "scale=3; (${homozygous_variant}-${exp_homo_var})^2/${exp_homo_var}" | bc -l)

echo -e ${coords}"\t"${missing}"\t"${homozygous_reference}"\t"${heterozygous}"\t"${homozygous_variant}"\t"${exp_homo_ref}"\t"${exp_het}"\t"${exp_homo_var}"\t"${chi_hom_ref}"\t"${chi_het}"\t"${chi_hom_var}; 
done < genotype_matrix >> genotype_matrix.chisquare













R
library(tidyverse)
library(zoo)

data_sus <- read.table("chr5_XQTL_susceptible_diff0.8_biased_variants.frq", header=F)
data_sus$population <- "XQTL_susceptible"

data_res <- read.table("chr5_XQTL_resistant_diff0.8_biased_variants.frq", header=F)
data_res$population <- "XQTL_resistant"

data <- bind_rows(data_sus, data_res)

a<- ggplot(data, aes(V2, V8, col=population)) +
    geom_point(size=0.1) + 
    facet_grid(.~V1) +
    theme_bw() 


data <- bind_rows(data_sus, data_res)
data <- data %>% mutate(rolling_avg = rollmean(V8, k=50, fill=NA, align='right'))
ggplot() + geom_line(aes(data$V2, data$rolling_avg, col=data$population), size=1) + ylim(0,1) + geom_vline(xintercept=31521884)


R
library(tidyverse)
library(zoo)

data_sus <- read.table("chr1_XQTL_susceptible_diff0.8_biased_variants.frq", header=F)
data_sus$population <- "XQTL_susceptible"

data_res <- read.table("chr1_XQTL_resistant_diff0.8_biased_variants.frq", header=F)
data_res$population <- "XQTL_resistant"

data <- bind_rows(data_sus, data_res)
data <- data %>% mutate(rolling_avg = rollmean(V8, k=10, fill=NA, align='right'))

ggplot(data, aes(V2, V8, col=population)) +
    geom_point(size=0.5) + 
    geom_smooth() +
    facet_grid(.~V1) +
    theme_bw() 


data <- data %>% mutate(rolling_avg = rollmean(V8, k=100, fill=NA, align='right'))
ggplot(data, aes(V2, rolling_avg, col=population)) + geom_line() + ylim(0,1)






vcftools \
    --gzvcf hcontortus_chr5_Celeg_TT_arrow_pilon.raw.vcf.gz \
    --keep XQTL_resistant.samples.keep_list \
    --positions MHCO3_v_MHCO18_diff0.8_biased_variants.keep-positions \
    --min-alleles 2 \
    --max-alleles 2 \
    --site-pi \
    --out chr5_XQTL_resistant_diff0.8_biased_variants


vcftools \
    --gzvcf hcontortus_chr5_Celeg_TT_arrow_pilon.raw.vcf.gz \
    --keep XQTL_susceptible.samples.keep_list \
    --positions MHCO3_v_MHCO18_diff0.8_biased_variants.keep-positions \
    --min-alleles 2 \
    --max-alleles 2 \
    --site-pi \
    --out chr5_XQTL_susceptible_diff0.8_biased_variants



pixy --vcf hcontortus_chr5_Celeg_TT_arrow_pilon.raw.vcf.gz --sites_file MHCO3_v_MHCO18_diff0.8_biased_variants.keep-positions --stats pi --populations XQTLsamples.keep_list.pixypop2 --bypass_invariant_check yes --window_size 1

pixy --vcf hcontortus_chr5_Celeg_TT_arrow_pilon.raw.vcf.gz --sites_file MHCO3_v_MHCO18_diff0.8_biased_variants.keep-positions --stats dxy --populations XQTLsamples.keep_list.pixypop2 --bypass_invariant_check yes --window_size 1