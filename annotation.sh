#Make the tst.awk file
BEGIN { OFS="\t" }
NR==FNR {
    ++cnt[$1]
    beg[$1,cnt[$1]] = $2
    end[$1,cnt[$1]] = $3
    val[$1,cnt[$1]] = $0
    next
}
$1 in cnt {
    for (i=1; i<=cnt[$1]; i++) {
        if ( (beg[$1,i] <= $2) && ($2 <= end[$1,i]) ) {
            print $0, val[$1,i]
        }
    }
}

#Prepare a modified gff file for each reference genome
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

chmod +x mod-gff.sh
./mod-gff.sh

#Replace contig names with numbers in your filtered vcf file, and only print specific columns of your vcf file 
#!/bin/bash
for i in $(ls *.vcf.filtered | rev | cut -c 14- | rev | uniq)
do
sed 's/NODE_/contig_/g' ${i}.vcf.filtered > ${i}.vcf.filtered.mod
sed 's/_length_[0-9]*_cov_[0-9]*.[0-9]*//g' ${i}.vcf.filtered.mod > ${i}.vcf.filtered.mod2
cat ${i}.vcf.filtered.mod2 | awk '{print $1, $2, $5}' > ${i}.vcf.filtered.mod3
rm ${i}.vcf.filtered.mod
rm ${i}.vcf.filtered.mod2
done 

chmod +x mod-vcf.sh
./mod-vcf.sh

#Move all modified files to the same directory

chmod +x tst.awk

#!/bin/bash
for i in $(ls 0944*.vcf.filtered.mod3 | rev | cut -c 19- | rev | uniq)
do
awk -f tst.awk AR-0944-mod4.gff ${i}.vcf.filtered.mod3 > ${i}.vcf.filtered.annotated
done
for i in $(ls 0945*.vcf.filtered.mod3 | rev | cut -c 19- | rev | uniq)
do
awk -f tst.awk AR-0945-mod4.gff ${i}.vcf.filtered.mod3 > ${i}.vcf.filtered.annotated
done
for i in $(ls 0948*.vcf.filtered.mod3 | rev | cut -c 19- | rev | uniq)
do
awk -f tst.awk AR-0948-mod4.gff ${i}.vcf.filtered.mod3 > ${i}.vcf.filtered.annotated
done
for i in $(ls 953*.vcf.filtered.mod3 | rev | cut -c 19- | rev | uniq)
do
awk -f tst.awk AR-0953-mod4.gff ${i}.vcf.filtered.mod3 > ${i}.vcf.filtered.annotated
done
for i in $(ls 0957*.vcf.filtered.mod3 | rev | cut -c 19- | rev | uniq)
do
awk -f tst.awk AR-0957-mod4.gff ${i}.vcf.filtered.mod3 > ${i}.vcf.filtered.annotated
done
for i in $(ls G1*.vcf.filtered.mod3 | rev | cut -c 19- | rev | uniq)
do
awk -f tst.awk AR-0953-mod4.gff ${i}.vcf.filtered.mod3 > ${i}.vcf.filtered.annotated
done
for i in $(ls G2*.vcf.filtered.mod3 | rev | cut -c 19- | rev | uniq)
do
awk -f tst.awk AR-0944-mod4.gff ${i}.vcf.filtered.mod3 > ${i}.vcf.filtered.annotated
done
for i in $(ls G3*.vcf.filtered.mod3 | rev | cut -c 19- | rev | uniq)
do
awk -f tst.awk AR-0957-mod4.gff ${i}.vcf.filtered.mod3 > ${i}.vcf.filtered.annotated
done
for i in $(ls G4*.vcf.filtered.mod3 | rev | cut -c 19- | rev | uniq)
do
awk -f tst.awk AR-0948-mod4.gff ${i}.vcf.filtered.mod3 > ${i}.vcf.filtered.annotated
done
for i in $(ls G5*.vcf.filtered.mod3 | rev | cut -c 19- | rev | uniq)
do
awk -f tst.awk AR-0945-mod4.gff ${i}.vcf.filtered.mod3 > ${i}.vcf.filtered.annotated
done

chmod +x annotation_fin.sh
./annotation_fin.sh

#!/bin/bash
for i in $(ls *.vcf.filtered.annotated | rev | cut -c 24- | rev | uniq)
do
sed 's/ /\t/g' ${i}.vcf.filtered.annotated > ${i}.vcf.filtered.annotated2
done

chmod +x whitespace.sh
./whitespace.sh

#Subtract mutations present in control replicated; write another file

#0944

C1 <- read.csv("0944-C1_S18_L001.vcf.filtered.annotated2", sep="\t", head=F)
C2 <- read.csv("0944-C2_S19_L001.vcf.filtered.annotated2", sep="\t", head=F)
C3 <- read.csv("0944-C3_S20_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen1 <- read.csv("0944-PEN-A_S14_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen2 <- read.csv("0944-PEN-B_S15_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen3 <- read.csv("0944-PEN-C_S16_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen4 <- read.csv("0944-PEN-D_S17_L001.vcf.filtered.annotated2", sep="\t", head=F)
Azi1 <- read.csv("G2S21_S5_L001.vcf.filtered.annotated2", sep="\t", head=F)
Azi2 <- read.csv("G2S41_S6_L001.vcf.filtered.annotated2", sep="\t", head=F)

C1 <- C1[,-4:-6]
C2 <- C2[,-4:-6]
C3 <- C3[,-4:-6]
Pen1 <- Pen1[,-4:-6]
Pen2 <- Pen2[,-4:-6]
Pen3 <- Pen3[,-4:-6]
Pen4 <- Pen4[,-4:-6]
Azi1 <- Azi1[,-4:-6]
Azi2 <- Azi2[,-4:-6]

C1$node_position <- with(C1, paste(V1, V2, sep="_"))
C2$node_position <- with(C2, paste(V1, V2, sep="_"))
C3$node_position <- with(C3, paste(V1, V2, sep="_"))
Pen1$node_position <- with(Pen1, paste(V1, V2, sep="_"))
Pen2$node_position <- with(Pen2, paste(V1, V2, sep="_"))
Pen3$node_position <- with(Pen3, paste(V1, V2, sep="_"))
Pen4$node_position <- with(Pen4, paste(V1, V2, sep="_"))
Azi1$node_position <- with(Azi1, paste(V1, V2, sep="_"))
Azi2$node_position <- with(Azi2, paste(V1, V2, sep="_"))

Incorrect <- rbind(C1,C2, C2)
Incorrect$x <- rep("x",dim(Incorrect)[1])
Incorrect <- Incorrect[,-1:-2]
wordList <- c("x")

Pen1 <- merge(Pen1, Incorrect, by="node_position", all.x=T, all.y=T)
Pen2 <- merge(Pen2, Incorrect, by="node_position", all.x=T, all.y=T)
Pen3 <- merge(Pen3, Incorrect, by="node_position", all.x=T, all.y=T)
Pen4 <- merge(Pen4, Incorrect, by="node_position", all.x=T, all.y=T)
Azi1 <- merge(Azi1, Incorrect, by="node_position", all.x=T, all.y=T)
Azi2 <- merge(Azi2, Incorrect, by="node_position", all.x=T, all.y=T)

Pen1 <- subset(Pen1, !(x %in% wordList))
Pen1$strain <- rep("Pen1", dim(Pen1)[1])
Pen2 <- subset(Pen2, !(x %in% wordList))
Pen2$strain <- rep("Pen2", dim(Pen2)[1])
Pen3 <- subset(Pen3, !(x %in% wordList))
Pen3$strain <- rep("Pen3", dim(Pen3)[1])
Pen4 <- subset(Pen4, !(x %in% wordList))
Pen4$strain <- rep("Pen4", dim(Pen4)[1])
Azi1 <- subset(Azi1, !(x %in% wordList))
Azi1$strain <- rep("Azi1", dim(Azi1)[1])
Azi2 <- subset(Azi2, !(x %in% wordList))
Azi2$strain <- rep("Azi2", dim(Azi2)[1])

write.table(Pen1, "~/Desktop/mut/0944_Pen1.txt")
write.table(Pen2, "~/Desktop/mut/0944_Pen2.txt")
write.table(Pen3, "~/Desktop/mut/0944_Pen3.txt")
write.table(Pen4, "~/Desktop/mut/0944_Pen4.txt")
write.table(Azi1, "~/Desktop/mut/0944_Azi1.txt")
write.table(Azi2, "~/Desktop/mut/0944_Azi2.txt")

#####
#Elongata

C1 <- read.csv("0945-C1_S15_L001.vcf.filtered.annotated2", sep="\t", head=F)
C2 <- read.csv("0945-C2_S16_L001.vcf.filtered.annotated2", sep="\t", head=F)
C3 <- read.csv("0945-C3_S17_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen1 <- read.csv("G5S11_S14_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen2 <- read.csv("G5S41_S16_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen3 <- read.csv("G5S31_S15_L001.vcf.filtered.annotated2", sep="\t", head=F)
#Pen4 <- read.csv("", sep="\t", head=F)
#Azi1 <- read.csv("", sep="\t", head=F)
Azi2 <- read.csv("0945-AZI-B_S18_L001.vcf.filtered.annotated2", sep="\t", head=F)
Azi3 <- read.csv("0945-AZI-C_S19_L001.vcf.filtered.annotated2", sep="\t", head=F)
Azi4 <- read.csv("0945-AZI-D_S20_L001.vcf.filtered.annotated2", sep="\t", head=F)

C1 <- C1[,-4:-6]
C2 <- C2[,-4:-6]
C3 <- C3[,-4:-6]
Pen1 <- Pen1[,-4:-6]
Pen2 <- Pen2[,-4:-6]
Pen3 <- Pen3[,-4:-6]
#Pen4 <- Pen4[,-4:-6]
#Azi1 <- Azi1[,-4:-6]
Azi2 <- Azi2[,-4:-6]
Azi3 <- Azi3[,-4:-6]
Azi4 <- Azi4[,-4:-6]

C1$node_position <- with(C1, paste(V1, V2, sep="_"))
C2$node_position <- with(C2, paste(V1, V2, sep="_"))
C3$node_position <- with(C3, paste(V1, V2, sep="_"))
Pen1$node_position <- with(Pen1, paste(V1, V2, sep="_"))
Pen2$node_position <- with(Pen2, paste(V1, V2, sep="_"))
Pen3$node_position <- with(Pen3, paste(V1, V2, sep="_"))
#Pen4$node_position <- with(Pen4, paste(V1, V2, sep="_"))
#Azi1$node_position <- with(Azi1, paste(V1, V2, sep="_"))
Azi2$node_position <- with(Azi2, paste(V1, V2, sep="_"))
Azi3$node_position <- with(Azi3, paste(V1, V2, sep="_"))
Azi4$node_position <- with(Azi4, paste(V1, V2, sep="_"))

Incorrect <- rbind(C1,C2, C2)
Incorrect$x <- rep("x",dim(Incorrect)[1])
Incorrect <- Incorrect[,-1:-2]
wordList <- c("x")

Pen1 <- merge(Pen1, Incorrect, by="node_position", all.x=T, all.y=T)
Pen2 <- merge(Pen2, Incorrect, by="node_position", all.x=T, all.y=T)
Pen3 <- merge(Pen3, Incorrect, by="node_position", all.x=T, all.y=T)
#Pen4 <- merge(Pen4, Incorrect, by="node_position", all.x=T, all.y=T)
#Azi1 <- merge(Azi1, Incorrect, by="node_position", all.x=T, all.y=T)
Azi2 <- merge(Azi2, Incorrect, by="node_position", all.x=T, all.y=T)
Azi3 <- merge(Azi3, Incorrect, by="node_position", all.x=T, all.y=T)
Azi4 <- merge(Azi4, Incorrect, by="node_position", all.x=T, all.y=T)

Pen1 <- subset(Pen1, !(x %in% wordList))
Pen1$strain <- rep("Pen1", dim(Pen1)[1])
Pen2 <- subset(Pen2, !(x %in% wordList))
Pen2$strain <- rep("Pen2", dim(Pen2)[1])
Pen3 <- subset(Pen3, !(x %in% wordList))
Pen3$strain <- rep("Pen3", dim(Pen3)[1])
#Pen4 <- subset(Pen4, !(x %in% wordList))
#Pen4$strain <- rep("Pen4", dim(Pen4)[1])
#Azi1 <- subset(Azi1, !(x %in% wordList))
#Azi1$strain <- rep("Azi1", dim(Azi1)[1])
Azi2 <- subset(Azi2, !(x %in% wordList))
Azi2$strain <- rep("Azi2", dim(Azi2)[1])
Azi3 <- subset(Azi3, !(x %in% wordList))
Azi3$strain <- rep("Azi3", dim(Azi2)[1])
Azi4 <- subset(Azi4, !(x %in% wordList))
Azi4$strain <- rep("Azi4", dim(Azi2)[1])

write.table(Pen1, "~/Desktop/mut/0945_Pen1.txt")
write.table(Pen2, "~/Desktop/mut/0945_Pen2.txt")
write.table(Pen3, "~/Desktop/mut/0945_Pen3.txt")
#write.csv(Pen4, "~/Desktop/mut/0945_Pen4.csv")
#write.csv(Azi1, "~/Desktop/mut/0945_Azi1.csv")
write.table(Azi2, "~/Desktop/mut/0945_Azi2.txt")
write.table(Azi3, "~/Desktop/mut/0945_Azi3.txt")
write.table(Azi4, "~/Desktop/mut/0945_Azi4.txt")

#####
#0948

C1 <- read.csv("0948-C1_S1_L001.vcf.filtered.annotated2", sep="\t", head=F)
C2 <- read.csv("0948-C2_S2_L001.vcf.filtered.annotated2", sep="\t", head=F)
C3 <- read.csv("0948-C3_S3_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen1 <- read.csv("G4S11_S11_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen2 <- read.csv("G4S12_S12_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen3 <- read.csv("G4S32_S13_L001.vcf.filtered.annotated2", sep="\t", head=F)
#Pen4 <- read.csv("", sep="\t", head=F)
Azi1 <- read.csv("0948-AZI-A_S10_L001.vcf.filtered.annotated2", sep="\t", head=F)
Azi2 <- read.csv("0948-AZI-B_S11_L001.vcf.filtered.annotated2", sep="\t", head=F)
Azi3 <- read.csv("0948-AZI-C_S12_L001.vcf.filtered.annotated2", sep="\t", head=F)
Azi4 <- read.csv("0948-AZI-D_S13_L001.vcf.filtered.annotated2", sep="\t", head=F)

C1 <- C1[,-4:-6]
C2 <- C2[,-4:-6]
C3 <- C3[,-4:-6]
Pen1 <- Pen1[,-4:-6]
Pen2 <- Pen2[,-4:-6]
Pen3 <- Pen3[,-4:-6]
#Pen4 <- Pen4[,-4:-6]
Azi1 <- Azi1[,-4:-6]
Azi2 <- Azi2[,-4:-6]
Azi3 <- Azi3[,-4:-6]
Azi4 <- Azi4[,-4:-6]

C1$node_position <- with(C1, paste(V1, V2, sep="_"))
C2$node_position <- with(C2, paste(V1, V2, sep="_"))
C3$node_position <- with(C3, paste(V1, V2, sep="_"))
Pen1$node_position <- with(Pen1, paste(V1, V2, sep="_"))
Pen2$node_position <- with(Pen2, paste(V1, V2, sep="_"))
Pen3$node_position <- with(Pen3, paste(V1, V2, sep="_"))
#Pen4$node_position <- with(Pen4, paste(V1, V2, sep="_"))
Azi1$node_position <- with(Azi1, paste(V1, V2, sep="_"))
Azi2$node_position <- with(Azi2, paste(V1, V2, sep="_"))
Azi3$node_position <- with(Azi3, paste(V1, V2, sep="_"))
Azi4$node_position <- with(Azi4, paste(V1, V2, sep="_"))

Incorrect <- rbind(C1,C2, C2)
Incorrect$x <- rep("x",dim(Incorrect)[1])
Incorrect <- Incorrect[,-1:-2]
wordList <- c("x")

Pen1 <- merge(Pen1, Incorrect, by="node_position", all.x=T, all.y=T)
Pen2 <- merge(Pen2, Incorrect, by="node_position", all.x=T, all.y=T)
Pen3 <- merge(Pen3, Incorrect, by="node_position", all.x=T, all.y=T)
#Pen4 <- merge(Pen4, Incorrect, by="node_position", all.x=T, all.y=T)
Azi1 <- merge(Azi1, Incorrect, by="node_position", all.x=T, all.y=T)
Azi2 <- merge(Azi2, Incorrect, by="node_position", all.x=T, all.y=T)
Azi3 <- merge(Azi3, Incorrect, by="node_position", all.x=T, all.y=T)
Azi4 <- merge(Azi4, Incorrect, by="node_position", all.x=T, all.y=T)

Pen1 <- subset(Pen1, !(x %in% wordList))
Pen1$strain <- rep("Pen1", dim(Pen1)[1])
Pen2 <- subset(Pen2, !(x %in% wordList))
Pen2$strain <- rep("Pen2", dim(Pen2)[1])
Pen3 <- subset(Pen3, !(x %in% wordList))
Pen3$strain <- rep("Pen3", dim(Pen3)[1])
#Pen4 <- subset(Pen4, !(x %in% wordList))
#Pen4$strain <- rep("Pen4", dim(Pen4)[1])
Azi1 <- subset(Azi1, !(x %in% wordList))
Azi1$strain <- rep("Azi1", dim(Azi1)[1])
Azi2 <- subset(Azi2, !(x %in% wordList))
Azi2$strain <- rep("Azi2", dim(Azi2)[1])
Azi3 <- subset(Azi3, !(x %in% wordList))
Azi3$strain <- rep("Azi3", dim(Azi2)[1])
Azi4 <- subset(Azi4, !(x %in% wordList))
Azi4$strain <- rep("Azi4", dim(Azi2)[1])

write.table(Pen1, "~/Desktop/mut/0948_Pen1.txt")
write.table(Pen2, "~/Desktop/mut/0948_Pen2.txt")
write.table(Pen3, "~/Desktop/mut/0948_Pen3.txt")
#write.csv(Pen4, "~/Desktop/mut/0948_Pen4.csv")
write.table(Azi1, "~/Desktop/mut/0948_Azi1.txt")
write.table(Azi2, "~/Desktop/mut/0948_Azi2.txt")
write.table(Azi3, "~/Desktop/mut/0948_Azi3.txt")
write.table(Azi4, "~/Desktop/mut/0948_Azi4.txt")

#####
#0953

C1 <- read.csv("953-C1_S7_L001.vcf.filtered.annotated2", sep="\t", head=F)
C2 <- read.csv("953-C2_S8_L001.vcf.filtered.annotated2", sep="\t", head=F)
C3 <- read.csv("953-C3_S9_L001.vcf.filtered.annotated2", sep="\t", head=F)
#Pen1 <- read.csv("G4S11_S11_L001.vcf.filtered.annotated2", sep="\t", head=F)
#Pen2 <- read.csv("G4S12_S12_L001.vcf.filtered.annotated2", sep="\t", head=F)
#Pen3 <- read.csv("G4S32_S13_L001.vcf.filtered.annotated2", sep="\t", head=F)
#Pen4 <- read.csv("", sep="\t", head=F)
Azi1 <- read.csv("G1S11_S1_L001.vcf.filtered.annotated2", sep="\t", head=F)
Azi2 <- read.csv("G1S21_S2_L001.vcf.filtered.annotated2", sep="\t", head=F)
Azi3 <- read.csv("G1S31_S3_L001.vcf.filtered.annotated2", sep="\t", head=F)
Azi4 <- read.csv("G1S41_S4_L001.vcf.filtered.annotated2", sep="\t", head=F)

C1 <- C1[,-4:-6]
C2 <- C2[,-4:-6]
C3 <- C3[,-4:-6]
#Pen1 <- Pen1[,-4:-6]
#Pen2 <- Pen2[,-4:-6]
#Pen3 <- Pen3[,-4:-6]
#Pen4 <- Pen4[,-4:-6]
Azi1 <- Azi1[,-4:-6]
Azi2 <- Azi2[,-4:-6]
Azi3 <- Azi3[,-4:-6]
Azi4 <- Azi4[,-4:-6]

C1$node_position <- with(C1, paste(V1, V2, sep="_"))
C2$node_position <- with(C2, paste(V1, V2, sep="_"))
C3$node_position <- with(C3, paste(V1, V2, sep="_"))
#Pen1$node_position <- with(Pen1, paste(V1, V2, sep="_"))
#Pen2$node_position <- with(Pen2, paste(V1, V2, sep="_"))
#Pen3$node_position <- with(Pen3, paste(V1, V2, sep="_"))
#Pen4$node_position <- with(Pen4, paste(V1, V2, sep="_"))
Azi1$node_position <- with(Azi1, paste(V1, V2, sep="_"))
Azi2$node_position <- with(Azi2, paste(V1, V2, sep="_"))
Azi3$node_position <- with(Azi3, paste(V1, V2, sep="_"))
Azi4$node_position <- with(Azi4, paste(V1, V2, sep="_"))

Incorrect <- rbind(C1,C2, C2)
Incorrect$x <- rep("x",dim(Incorrect)[1])
Incorrect <- Incorrect[,-1:-2]
wordList <- c("x")

#Pen1 <- merge(Pen1, Incorrect, by="node_position", all.x=T, all.y=T)
#Pen2 <- merge(Pen2, Incorrect, by="node_position", all.x=T, all.y=T)
#Pen3 <- merge(Pen3, Incorrect, by="node_position", all.x=T, all.y=T)
#Pen4 <- merge(Pen4, Incorrect, by="node_position", all.x=T, all.y=T)
Azi1 <- merge(Azi1, Incorrect, by="node_position", all.x=T, all.y=T)
Azi2 <- merge(Azi2, Incorrect, by="node_position", all.x=T, all.y=T)
Azi3 <- merge(Azi3, Incorrect, by="node_position", all.x=T, all.y=T)
Azi4 <- merge(Azi4, Incorrect, by="node_position", all.x=T, all.y=T)

#Pen1 <- subset(Pen1, !(x %in% wordList))
#Pen1$strain <- rep("Pen1", dim(Pen1)[1])
#Pen2 <- subset(Pen2, !(x %in% wordList))
#Pen2$strain <- rep("Pen2", dim(Pen2)[1])
#Pen3 <- subset(Pen3, !(x %in% wordList))
#Pen3$strain <- rep("Pen3", dim(Pen3)[1])
#Pen4 <- subset(Pen4, !(x %in% wordList))
#Pen4$strain <- rep("Pen4", dim(Pen4)[1])
Azi1 <- subset(Azi1, !(x %in% wordList))
Azi1$strain <- rep("Azi1", dim(Azi1)[1])
Azi2 <- subset(Azi2, !(x %in% wordList))
Azi2$strain <- rep("Azi2", dim(Azi2)[1])
Azi3 <- subset(Azi3, !(x %in% wordList))
Azi3$strain <- rep("Azi3", dim(Azi2)[1])
Azi4 <- subset(Azi4, !(x %in% wordList))
Azi4$strain <- rep("Azi4", dim(Azi2)[1])

write.table(Azi1, "~/Desktop/mut/0953_Azi1.txt")
write.table(Azi2, "~/Desktop/mut/0953_Azi2.txt")
write.table(Azi3, "~/Desktop/mut/0953_Azi3.txt")
write.table(Azi4, "~/Desktop/mut/0953_Azi4.txt")

#####
#0957

C1 <- read.csv("0957-C1_S4_L001.vcf.filtered.annotated2", sep="\t", head=F)
C2 <- read.csv("0957-C2_S5_L001.vcf.filtered.annotated2", sep="\t", head=F)
C3 <- read.csv("0957-C3_S6_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen1 <- read.csv("G3S11_S7_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen2 <- read.csv("G3S21_S8_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen3 <- read.csv("G3S31_S9_L001.vcf.filtered.annotated2", sep="\t", head=F)
Pen4 <- read.csv("G3S41_S10_L001.vcf.filtered.annotated2", sep="\t", head=F)
#Azi1 <- read.csv("G1S11_S1_L001.vcf.filtered.annotated2", sep="\t", head=F)
#Azi2 <- read.csv("G1S21_S2_L001.vcf.filtered.annotated2", sep="\t", head=F)
#Azi3 <- read.csv("G1S31_S3_L001.vcf.filtered.annotated2", sep="\t", head=F)
#Azi4 <- read.csv("G1S41_S4_L001.vcf.filtered.annotated2", sep="\t", head=F)

C1 <- C1[,-4:-6]
C2 <- C2[,-4:-6]
C3 <- C3[,-4:-6]
Pen1 <- Pen1[,-4:-6]
Pen2 <- Pen2[,-4:-6]
Pen3 <- Pen3[,-4:-6]
Pen4 <- Pen4[,-4:-6]
#Azi1 <- Azi1[,-4:-6]
#Azi2 <- Azi2[,-4:-6]
#Azi3 <- Azi3[,-4:-6]
#Azi4 <- Azi4[,-4:-6]

C1$node_position <- with(C1, paste(V1, V2, sep="_"))
C2$node_position <- with(C2, paste(V1, V2, sep="_"))
C3$node_position <- with(C3, paste(V1, V2, sep="_"))
Pen1$node_position <- with(Pen1, paste(V1, V2, sep="_"))
Pen2$node_position <- with(Pen2, paste(V1, V2, sep="_"))
Pen3$node_position <- with(Pen3, paste(V1, V2, sep="_"))
Pen4$node_position <- with(Pen4, paste(V1, V2, sep="_"))
#Azi1$node_position <- with(Azi1, paste(V1, V2, sep="_"))
#Azi2$node_position <- with(Azi2, paste(V1, V2, sep="_"))
#Azi3$node_position <- with(Azi3, paste(V1, V2, sep="_"))
#Azi4$node_position <- with(Azi4, paste(V1, V2, sep="_"))

Incorrect <- rbind(C1,C2, C2)
Incorrect$x <- rep("x",dim(Incorrect)[1])
Incorrect <- Incorrect[,-1:-2]
wordList <- c("x")

Pen1 <- merge(Pen1, Incorrect, by="node_position", all.x=T, all.y=T)
Pen2 <- merge(Pen2, Incorrect, by="node_position", all.x=T, all.y=T)
Pen3 <- merge(Pen3, Incorrect, by="node_position", all.x=T, all.y=T)
Pen4 <- merge(Pen4, Incorrect, by="node_position", all.x=T, all.y=T)
#Azi1 <- merge(Azi1, Incorrect, by="node_position", all.x=T, all.y=T)
#Azi2 <- merge(Azi2, Incorrect, by="node_position", all.x=T, all.y=T)
#Azi3 <- merge(Azi3, Incorrect, by="node_position", all.x=T, all.y=T)
#Azi4 <- merge(Azi4, Incorrect, by="node_position", all.x=T, all.y=T)

Pen1 <- subset(Pen1, !(x %in% wordList))
Pen1$strain <- rep("Pen1", dim(Pen1)[1])
Pen2 <- subset(Pen2, !(x %in% wordList))
Pen2$strain <- rep("Pen2", dim(Pen2)[1])
Pen3 <- subset(Pen3, !(x %in% wordList))
Pen3$strain <- rep("Pen3", dim(Pen3)[1])
Pen4 <- subset(Pen4, !(x %in% wordList))
Pen4$strain <- rep("Pen4", dim(Pen4)[1])
#Azi1 <- subset(Azi1, !(x %in% wordList))
#Azi1$strain <- rep("Azi1", dim(Azi1)[1])
#Azi2 <- subset(Azi2, !(x %in% wordList))
#Azi2$strain <- rep("Azi2", dim(Azi2)[1])
#Azi3 <- subset(Azi3, !(x %in% wordList))
#Azi3$strain <- rep("Azi3", dim(Azi2)[1])
#Azi4 <- subset(Azi4, !(x %in% wordList))
#Azi4$strain <- rep("Azi4", dim(Azi2)[1])

write.table(Pen1, "~/Desktop/mut/0957_Pen1.txt")
write.table(Pen2, "~/Desktop/mut/0957_Pen2.txt")
write.table(Pen3, "~/Desktop/mut/0957_Pen3.txt")
write.table(Pen4, "~/Desktop/mut/0957_Pen4.txt")

sed 's/"//g' 0944_Azi1.txt > 0944_Azi1.txt2
sed 's/"//g' 0944_Azi1.txt > 0944_Azi1.txt2
sed 's/"//g' 0944_Pen2.txt > 0944_Pen2.txt2
sed 's/"//g' 0945_Azi2.txt > 0945_Azi2.txt2
sed 's/"//g' 0945_Pen1.txt > 0945_Pen1.txt2
sed 's/"//g' 0948_Azi1.txt > 0948_Azi1.txt2
sed 's/"//g' 0948_Azi4.txt > 0948_Azi4.txt2
sed 's/"//g' 0948_Pen3.txt > 0948_Pen3.txt2
sed 's/"//g' 0953_Azi3.txt > 0953_Azi3.txt2
sed 's/"//g' 0957_Pen2.txt > 0957_Pen2.txt2
sed 's/"//g' 0944_Azi2.txt > 0944_Azi2.txt2
sed 's/"//g' 0944_Pen3.txt > 0944_Pen3.txt2
sed 's/"//g' 0945_Azi3.txt > 0945_Azi3.txt2
sed 's/"//g' 0945_Pen2.txt > 0945_Pen2.txt2
sed 's/"//g' 0948_Azi2.txt > 0948_Azi2.txt2
sed 's/"//g' 0948_Pen1.txt > 0948_Pen1.txt2
sed 's/"//g' 0953_Azi1.txt > 0953_Azi1.txt2
sed 's/"//g' 0953_Azi4.txt > 0953_Azi4.txt2
sed 's/"//g' 0957_Pen3.txt > 0957_Pen3.txt2
sed 's/"//g' 0944_Pen1.txt > 0944_Pen1.txt2
sed 's/"//g' 0944_Pen4.txt > 0944_Pen4.txt2
sed 's/"//g' 0945_Azi4.txt > 0945_Azi4.txt2
sed 's/"//g' 0945_Pen3.txt > 0945_Pen3.txt2
sed 's/"//g' 0948_Azi3.txt > 0948_Azi3.txt2
sed 's/"//g' 0948_Pen2.txt > 0948_Pen2.txt2
sed 's/"//g' 0953_Azi2.txt > 0953_Azi2.txt2
sed 's/"//g' 0957_Pen1.txt > 0957_Pen1.txt2
sed 's/"//g' 0957_Pen4.txt > 0957_Pen4.txt2

#Grab product information only, use BBedit
Find:ID=.+product=(.+) NA NA NA 
Replace:\1 

Find:;Dbxref.+ 
Replace: 

Find:node_position V1 V2 V3.x V7.x V3.y V7.y x strain
Replace:line node_position contig position mutation product drug

#Within lists concatenate to only one occurrence
#By drug, pen


A1 <- read.csv("0944_Pen1.txt2",sep="\t", head=T)
A2<- read.csv("0944_Pen2.txt2",sep="\t", head=T)
A3<- read.csv("0944_Pen3.txt2",sep="\t", head=T)
A4<- read.csv("0944_Pen4.txt2",sep="\t", head=T)

list1 <- c(unique(A1$product), unique(A2$product), unique(A3$product), unique(A4$product))
df1 <- data.frame(list1, rep("0944", length(list1)))
colnames(df1) <- c("product", "species")

B1<- read.csv("0945_Pen1.txt2",sep="\t", head=T)
B2<- read.csv("0945_Pen2.txt2",sep="\t", head=T)
B3<- read.csv("0945_Pen3.txt2",sep="\t", head=T)

list2 <-c(unique(B1$product), unique(B2$product), unique(B3$product))
df2 <- data.frame(list2, rep("0945", length(list2)))
colnames(df2) <- c("product", "species")

C1<- read.csv("0948_Pen1.txt2",sep="\t", head=T)
C2<- read.csv("0948_Pen2.txt2",sep="\t", head=T)
C3<- read.csv("0948_Pen3.txt2",sep="\t", head=T)

list3 <-c(unique(C1$product), unique(C2$product), unique(C3$product))
df3 <- data.frame(list3, rep("0948", length(list3)))
colnames(df3) <- c("product", "species")

D1<- read.csv("0957_Pen1.txt2",sep="\t", head=T)
D2<- read.csv("0957_Pen2.txt2",sep="\t", head=T)
D3<- read.csv("0957_Pen3.txt2",sep="\t", head=T)
D4<- read.csv("0957_Pen4.txt2",sep="\t", head=T)

list4 <-c(unique(D1$product), unique(D2$product), unique(D3$product), unique(D4$product))
df4 <- data.frame(list4, rep("0957", length(list4)))
colnames(df4) <- c("product", "species")

library(dplyr)
df1b <- df1 %>%
  group_by(product) %>%
  add_count(name = "count")
df2b <- df2 %>%
  group_by(product) %>%
  add_count(name = "count")
df3b <- df3 %>%
  group_by(product) %>%
  add_count(name = "count")
df4b <- df4 %>%
  group_by(product) %>%
  add_count(name = "count")

df1c <- df1b %>% distinct()
df2c <- df2b %>% distinct()
df3c <- df3b %>% distinct()
df4c <- df4b %>% distinct()

fin <- rbind(df1c,df2c,df3c,df4c)
ggplot(fin, aes(product, species)) + geom_tile(aes(fill = count))

#Modify table to contain 0 counts
write.csv(fin, "~/Desktop/fin_pen.csv")
fin2 <- read.csv("~/Desktop/fin_pen.csv")
ggplot(fin2, aes(product, as.character(species))) + geom_tile(aes(fill = count)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +theme(axis.text.x = element_text(face="italic"))

#By drug, azi

A1 <- read.csv("0944_Azi1.txt2",sep="\t", head=T)
A2<- read.csv("0944_Azi2.txt2",sep="\t", head=T)

list1 <- c(unique(A1$product), unique(A2$product))
df1 <- data.frame(list1, rep("0944", length(list1)))
colnames(df1) <- c("product", "species")

B2<- read.csv("0945_Azi2.txt2",sep="\t", head=T)
B3<- read.csv("0945_Azi3.txt2",sep="\t", head=T)
B4<- read.csv("0945_Azi4.txt2",sep="\t", head=T)

list2 <-c(unique(B2$product), unique(B3$product), unique(B4$product))
df2 <- data.frame(list2, rep("0945", length(list2)))
colnames(df2) <- c("product", "species")

C1<- read.csv("0948_Azi1.txt2",sep="\t", head=T)
C2<- read.csv("0948_Azi2.txt2",sep="\t", head=T)
C3<- read.csv("0948_Azi3.txt2",sep="\t", head=T)
C4<- read.csv("0948_Azi4.txt2",sep="\t", head=T)

list3 <-c(unique(C1$product), unique(C2$product), unique(C3$product), unique(C4$product))
df3 <- data.frame(list3, rep("0948", length(list3)))
colnames(df3) <- c("product", "species")

D1<- read.csv("0953_Azi1.txt2",sep="\t", head=T)
D2<- read.csv("0953_Azi2.txt2",sep="\t", head=T)
D3<- read.csv("0953_Azi3.txt2",sep="\t", head=T)
D4<- read.csv("0953_Azi4.txt2",sep="\t", head=T)

list4 <-c(unique(D1$product), unique(D2$product), unique(D3$product), unique(D4$product))
df4 <- data.frame(list4, rep("0953", length(list4)))
colnames(df4) <- c("product", "species")

library(dplyr)
df1b <- df1 %>%
  group_by(product) %>%
  add_count(name = "count")
df2b <- df2 %>%
  group_by(product) %>%
  add_count(name = "count")
df3b <- df3 %>%
  group_by(product) %>%
  add_count(name = "count")
df4b <- df4 %>%
  group_by(product) %>%
  add_count(name = "count")

df1c <- df1b %>% distinct()
df2c <- df2b %>% distinct()
df3c <- df3b %>% distinct()
df4c <- df4b %>% distinct()

fin <- rbind(df1c,df2c,df3c,df4c)
ggplot(fin, aes(product, species)) + geom_tile(aes(fill = count))

#Modify table to contain 0 counts
write.csv(fin, "~/Desktop/fin_azi.csv")
fin2 <- read.csv("~/Desktop/fin_azi.csv")
ggplot(fin2, aes(product, as.character(species))) + geom_tile(aes(fill = count)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +theme(axis.text.x = element_text(face="italic"))




