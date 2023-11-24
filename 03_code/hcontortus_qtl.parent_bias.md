# H.contortus QTL: Parental bias 



```bash 

cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/QTL/05_ANALYSIS/PARENT_BIAS


ln -s ../../QTL_MAPPING/MHCO3_P0_L3_n200_01/MHCO3_P0_L3_n200_01.bam
ln -s ../../QTL_MAPPING/MHCO3_P0_L3_n200_01/MHCO3_P0_L3_n200_01.bam.bai
ln -s ../../QTL_MAPPING/MHCO18_P0_L3_n200_IVM_01/MHCO18_P0_L3_n200_IVM_01.bam
ln -s ../../QTL_MAPPING/MHCO18_P0_L3_n200_IVM_01/MHCO18_P0_L3_n200_IVM_01.bam.bai

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

```


```bash
# median coverage of each group
head -n 10000000 MHCO3_v_MHCO18_poolsfrequency.csv | datamash --header-in median 5,7

#> 58	72
```


```bash
# checking variants with biased allele frequency
#-- freq diff = 0.8
#-- min coverage = 0.33 * median
awk '{print $1,$2,$8,$7,$6,$5,$6-$8}' OFS="\t" MHCO3_v_MHCO18_poolsfrequency.csv |\
    awk '{if($7>0.8 && $4>23 && $6>19) print}' OFS="\t" |\
    cut -f1 | sort | uniq -c


# 0.8 cutoff - 35711 total
   7024 hcontortus_chr1_Celeg_TT_arrow_pilon
  12695 hcontortus_chr2_Celeg_TT_arrow_pilon
   2247 hcontortus_chr3_Celeg_TT_arrow_pilon
   4808 hcontortus_chr4_Celeg_TT_arrow_pilon
   8244 hcontortus_chr5_Celeg_TT_arrow_pilon
    693 hcontortus_chrX_Celeg_TT_arrow_pilon 

# also tested a 0.9 cutoff, to check variant numbers
# 0.9 cutoff - 10706 total
   2166 hcontortus_chr1_Celeg_TT_arrow_pilon
   3721 hcontortus_chr2_Celeg_TT_arrow_pilon
    568 hcontortus_chr3_Celeg_TT_arrow_pilon
   1282 hcontortus_chr4_Celeg_TT_arrow_pilon
   2811 hcontortus_chr5_Celeg_TT_arrow_pilon
    158 hcontortus_chrX_Celeg_TT_arrow_pilon

#CHR POS MHCO3_freq MHCO3_cov MHCO18_freq MHCO18_cov diff(MHCO18_freq-MHCO3_freq)
```


## checking 
```bash
awk '{print $1,$2,$8,$7,$6,$5,$6-$8}' OFS="\t" MHCO3_v_MHCO18_poolsfrequency.csv |    awk '{if($7>0.8 && $4>23 && $6>19) print $1,$2 $7}' OFS="\t" > MHCO3_v_MHCO18_poolsfrequency.0.8freq.pos
```

```R
library(tidyverse)

data <- read.table("MHCO3_v_MHCO18_poolsfrequency.0.8freq.pos", header=F)

ggplot(data, aes(V2*1e6,V3,col=V1)) + 
    geom_point(size=0.1) + 
    facet_grid(V1~.) + 
    ylim(0,1) +
    theme_bw() +
    guides(color = FALSE) +
    labs(title="Variant positions that differ between MHco3 and MHco18 by > 0.8", 
        y="Variant frequency difference (min: MHCO18_freq-MHCO3_freq >0.8)", 
        x="Chromosomal position")

ggsave("MHCO3_v_MHCO18_poolsfrequency.0.8freq.pos.pdf", height=7, width=7)
ggsave("MHCO3_v_MHCO18_poolsfrequency.0.8freq.pos.png")
```
![](../04_analysis/MHCO3_v_MHCO18_poolsfrequency.0.8freq.pos.png)