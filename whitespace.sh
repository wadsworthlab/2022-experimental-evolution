#!/bin/bash
for i in $(ls *.vcf.filtered.annotated | rev | cut -c 24- | rev | uniq)
do
sed 's/ /\t/g' ${i}.vcf.filtered.annotated > ${i}.vcf.filtered.annotated2
done
