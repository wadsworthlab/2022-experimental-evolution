#!/bin/bash
for i in $(ls *.vcf.filtered | rev | cut -c 14- | rev | uniq)
do
sed 's/NODE_/contig_/g' ${i}.vcf.filtered > ${i}.vcf.filtered.mod
sed 's/_length_[0-9]*_cov_[0-9]*.[0-9]*//g' ${i}.vcf.filtered.mod > ${i}.vcf.filtered.mod2
cat ${i}.vcf.filtered.mod2 | awk '{print $1, $2, $5}' > ${i}.vcf.filtered.mod3
rm ${i}.vcf.filtered.mod
rm ${i}.vcf.filtered.mod2
done 
