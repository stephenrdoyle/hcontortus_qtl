


## plink

```bash

module load plink/1.90b6.18--h516909a_0

cd /lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL/05_ANALYSIS/PLINK


# use plink to import the vcf - it will write outputs
plink --vcf hcontortus_chr5_Celeg_TT_arrow_pilon.raw.vcf.gz --double-id --allow-extra-chr

# once it has written tandard outputs, use "--bfile" to input
plink --bfile plink --hardy --allow-extra-chr

--allow-no-sex 

--pheno file
```