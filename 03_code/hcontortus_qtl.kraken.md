



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