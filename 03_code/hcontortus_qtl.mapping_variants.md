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








## Variant calling



```bash
cd /lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL

find ~+ -type f -name '*.bam' | sort -V | grep -v "BLANK" | grep -v "Control" > 04_VARIANTS/bam.list

cd /lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL/04_VARIANTS/

./run_gatk.sh
```














#run_genotyping.sh 

export PREFIX=HCON_QTL  # prefix for output files
export REFERENCE=/nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/QTL/01_REFERENCE/HAEM_V4_final.chr.fa  # path to reference genome
export BAM_LIST=/lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL/04_VARIANTS/bam.list2  # path to list of BAM files

# Load GATK module
module load gatk/4.1.4.1


export WD=/lustre/scratch125/pam/teams/team333/sd21/haemonchus_contortus/QTL/04_VARIANTS/gatk_hc_HCON_QTL

# Define file locations
export LOG_FILES="${WD}/LOG_FILES"  # directory for log files
export REFERENCE_FILES="${WD}/REFERENCE_FILES"  # directory for reference files
export GATK_HC_GVCFs="${WD}/GATK_HC_GVCFs"  # directory for GATK HC GVCF files
export GATK_HC_MERGED="${WD}/GATK_HC_MERGED"  # directory for merged haplotype caller files

# Create directories if they don't exist
[ -d ${LOG_FILES} ] || mkdir -p ${LOG_FILES}
[ -d ${REFERENCE_FILES} ] || mkdir -p ${REFERENCE_FILES}
[ -d ${GATK_HC_GVCFs} ] || mkdir -p ${GATK_HC_GVCFs}
[ -d ${GATK_HC_MERGED} ] || mkdir -p ${GATK_HC_MERGED}



# Save current script in run folder to reproduce the exact output
cp ${PWD}/run_gatk_hc.sh ${PWD}/gatk_hc_${PREFIX}/commands.$(date -Iminutes).txt


#-------------------------------------------------------------------------------
### 03. GenomicsDBImport
#-------------------------------------------------------------------------------
func_GenomicsDBImport() {

ls -1 ${GATK_HC_GVCFs}/*complete/*gz > ${GATK_HC_MERGED}/gvcf.list

[ -d ${GATK_HC_MERGED}/LOGFILES ] || mkdir -p ${GATK_HC_MERGED}/LOGFILES


n=1
for SEQUENCE in ${REFERENCE_FILES}/REFsplit*list; do

    SEQUENCE=$( echo ${SEQUENCE} | awk -F '/' '{print $NF}' )

     mkdir ${GATK_HC_MERGED}/.tmp_${SEQUENCE}

    echo -e "gatk GenomicsDBImport --genomicsdb-workspace-path genomicsdb_${SEQUENCE} -R ${REFERENCE_FILES}/REF.fa --intervals ${REFERENCE_FILES}/${SEQUENCE} --reader-threads 20 --batch-size 100 \\" > ${GATK_HC_MERGED}/run_GenomicsDBImport.tmp.job_${n}
    while read SAMPLE; do
        echo -e "--variant ${SAMPLE} \\" >> ${GATK_HC_MERGED}/run_GenomicsDBImport.tmp.job_${n};
   done < ${GATK_HC_MERGED}/gvcf.list
   echo -e "--tmp-dir ${GATK_HC_MERGED}/.tmp_${SEQUENCE}" >> ${GATK_HC_MERGED}/run_GenomicsDBImport.tmp.job_${n};
   let "n+=1";
done

chmod a+x ${GATK_HC_MERGED}/run_GenomicsDBImport.tmp.job_*

# setup job conditions
JOBS=$( ls -1 ${GATK_HC_MERGED}/run_GenomicsDBImport.tmp.job_* | wc -l )
ID="U$(date +%s)"

#submit job array to run GenomicsDBImport
bsub -q long -R'span[hosts=1] select[mem>30000] rusage[mem=30000]' -n 20 -M30000 -J "gatk_GenomicsDBImport_[1-$JOBS]%100" -e "${GATK_HC_MERGED}/LOGFILES/gatk_GenomicsDBImport_[1-$JOBS].e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_GenomicsDBImport_[1-$JOBS].o" "${GATK_HC_MERGED}/run_GenomicsDBImport.tmp.job_\$LSB_JOBINDEX"

rm ${GATK_HC_MERGED}/MERGE_ARRAY_FINISHED
bsub -w "done(gatk_GenomicsDBImport_)" -q normal -R'span[hosts=1] select[mem>100] rusage[mem=100]' -n 1 -M100 -J "gatk_GenomicsDBImport_finish" -e "${GATK_HC_MERGED}/LOGFILES/gatk_GenomicsDBImport_finish.e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_GenomicsDBImport_finish.o" "touch ${GATK_HC_MERGED}/GenomicsDBImport_FINISHED"

until [ -f "${GATK_HC_MERGED}/GenomicsDBImport_FINISHED" ]
do
     sleep 10
done

}

export -f func_GenomicsDBImport

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
### 04. Genotype GVCFs
#-------------------------------------------------------------------------------

func_genotype_gvcfs() {

# split each chromosome up into separate jobs, and run genotyping on each individually.
n=1
for SEQUENCE in ${REFERENCE_FILES}/REFsplit*list; do
    SEQUENCE=$( echo ${SEQUENCE} | awk -F '/' '{print $NF}' )
    echo -e "gatk GenotypeGVCFs \
    -R ${REFERENCE_FILES}/REF.fa \
    -V gendb://genomicsdb_${SEQUENCE} \
    -O ${GATK_HC_MERGED}/${SEQUENCE}.cohort.tmp.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation" > ${GATK_HC_MERGED}/run_hc_genotype.tmp.job_${n};
    let "n+=1";
done

chmod a+x ${GATK_HC_MERGED}/run_hc_genotype*

# setup job conditions
JOBS=$( ls -1 ${GATK_HC_MERGED}/run_hc_genotype* | wc -l )
ID="U$(date +%s)"

bsub -q long -R'span[hosts=1] select[mem>20000] rusage[mem=20000]' -n 6 -M20000 -J "gatk_genotype_cohort_gvcf_[1-$JOBS]" -e "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_[1-$JOBS].e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_[1-$JOBS].o" "${GATK_HC_MERGED}/run_hc_genotype.tmp.job_*\$LSB_JOBINDEX"

rm ${GATK_HC_MERGED}/GENOTYPE_ARRAY_FINISHED
bsub -w "done(gatk_genotype_cohort_gvcf_)" -q normal -R'span[hosts=1] select[mem>100] rusage[mem=100]' -n 1 -M100 -J "gatk_genotype_cohort_gvcf_finish" -e "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_finish.e" -o "${GATK_HC_MERGED}/LOGFILES/gatk_genotype_cohort_gvcf_finish.o" "touch ${GATK_HC_MERGED}/GENOTYPE_ARRAY_FINISHED"

until [ -f "${GATK_HC_MERGED}/GENOTYPE_ARRAY_FINISHED" ]
do
     sleep 10
done

}

export -f func_genotype_gvcfs



#-------------------------------------------------------------------------------
### 05. Finish making VCF and cleanup
#-------------------------------------------------------------------------------


func_finish_vcf() {

    #
    ls ${GATK_HC_MERGED}/*.cohort.tmp.vcf.gz > ${GATK_HC_MERGED}/cohort.vcf.list

    # concatenate the vcf files in the list
    vcf-concat --files ${GATK_HC_MERGED}/cohort.vcf.list > ${GATK_HC_MERGED}/${PREFIX}.cohort.$(date -I).vcf

    # Compress the combined VCF file with bgzip
    bgzip -f ${GATK_HC_MERGED}/${PREFIX}.cohort.$(date -I).vcf

    # Create a tabix index for the compressed combined VCF file
    tabix -f ${GATK_HC_MERGED}/${PREFIX}.cohort.$(date -I).vcf.gz

    # Remove all files in the directory specified by GATK_HC_MERGED that match the pattern *tmp*
    rm ${GATK_HC_MERGED}/*tmp*

}

export -f func_finish_vcf





#-------------------------------------------------------------------------------
# running the pipeline
#-------------------------------------------------------------------------------




# func_merge_gvcf
bsub -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/gatk_03_GenomicsDBImport.o -e ${LOG_FILES}/gatk_03_GenomicsDBImport.e -J gatk_03_GenomicsDBImport_${PREFIX} func_GenomicsDBImport

# func_genotype_gvcfs
bsub -w "done(gatk_03_GenomicsDBImport_${PREFIX})" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>50000] rusage[mem=50000]" -q long -M50000 -n20 -o ${LOG_FILES}/gatk_04_genotype_gvcfs.o -e ${LOG_FILES}/gatk_04_genotype_gvcfs.e -J gatk_04_genotype_gvcfs_${PREFIX} func_genotype_gvcfs

# func_finish_vcf
bsub -w "done(gatk_04_genotype_gvcfs_${PREFIX})" -E 'test -e /nfs/users/nfs_s/sd21' -R "select[mem>1000] rusage[mem=1000]" -q long -M1000 -n1 -o ${LOG_FILES}/gatk_05_finish_vcf.o -e ${LOG_FILES}/gatk_05_finish_vcf.e -J gatk_05_finish_vcf_${PREFIX} func_finish_vcf

