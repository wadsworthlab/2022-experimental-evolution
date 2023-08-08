#!/bin/bash
for i in $(ls *.gff3 | rev | cut -c 6- | rev | uniq)
do
grep '^contig' ${i}.gff3 > ${i}-mod.gff
grep -v 'Bakta\tregion' ${i}-mod.gff > ${i}-mod2.gff
cat ${i}-mod2.gff | tr " " "-" > ${i}-mod3.gff
cat ${i}-mod3.gff | awk '{print $1, $4, $5, $9}' > ${i}-mod4.gff
rm  ${i}-mod.gff
rm ${i}-mod2.gff
rm ${i}-mod3.gff
done

