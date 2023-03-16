
input=${input:-"${set}/hum/${feat}.hg38"}
output=${output:-"${set}/dog/${feat}.hg38-to-cf3"}

${bedtools} sort -i ${input}.bed | \
${bedtools} merge -i - -c 4,5,6 -o distinct -d 10



cd ${dir}

# gff to bed
# /seq/vgb/software/bedops/bin/gff2bed < ${set}/hum/${set}.hg38.bed > ${set}/hum/${set}.hg38.sort.bed

# sort
# /seq/vgb

# re-map: hum to dog
# CrossMap.py bed ${hum2dog} ${set}/hum/${set}.hg38.sort.bed ${set}/dog/${set}.hg38-to-cf3.bed

# sort
 ${bedtools} sort -i ${set}/dog/${set}.hg38-to-cf3.bed > ${set}/dog/${set}.hg38-to-cf3.sort.bed

# merge
 ${bedtools} merge -i ${set}/dog/${set}.hg38-to-cf3.sort.bed -c 4,5,6 -o distinct -d 10 > ${set}/dog/${set}.hg38-to-cf3.sort.merge.bed

# back-map: dog to hum
# CrossMap.py bed ${dog2hum} ${set}/dog/${set}.hg38-to-cf3.sort.merge.bed ${set}/hum/${set}.hg38-to-cf3-to-hg38.bed

# sort
# ${bedtools} sort -i ${set}/hum/${set}.hg38-to-cf3-to-hg38.bed > ${set}/hum/${set}.hg38-to-cf3-to-hg38.sort.bed

# check reciprocal: intersect
${bedtools} intersect -wa -a ${set}/hum/${set}.hg38.sort.bed -b ${set}/hum/${set}.hg38-to-cf3-to-hg38.sort.bed | awk -F "\t" '{print $5}' | sort | uniq > ${set}/hum/${set}.hg38-to-cf3-to-hg38.recip.txt

# keep reciprocal
grep -wf ${set}/hum/${set}.hg38-to-cf3-to-hg38.recip.txt ${set}/dog/${set}.hg38-to-cf3.sort.merge.bed > ${set}/dog/${set}.hg38-to-cf3.sort.merge.recip.bed

# add flank and merge
${bedtools} sort -i ${set}/dog/${set}.hg38-to-cf3.sort.merge.recip.bed | ${bedtools} slop -g cf3.chrom.sizes -b 1000 -i - | ${bedtools} merge -i - -c 4,5,6 -o distinct -d 10 > ${set}/dog/${set}.hg38-to-cf3.sort.merge.recip.final.bed

# get matching
awk -F "\t" '{print $5}' ${set}/dog/${set}.hg38-to-cf3.shared.bed | sort | uniq > ${set}/dog/${set}.hg38-to-cf3.shared.txt

# final
grep -wf ${set}/dog/${set}.hg38-to-cf3.shared.txt ${set}/hum/${set}.hg38.sort.bed | ${bedtools} slop -g hg38.chrom.sizes -b 1000 -i - | ${bedtools} merge -i - -c 4,5,6 -o distinct -d 10 > ${set}/hum/${set}.hg38.shared.bed





# gff to bed
# /seq/vgb/software/bedops/bin/gff2bed < ${set}/hum/${feat}.hg38.gff3.gz > ${set}/hum/${feat}.hg38.split.bed

# sort
${bedtools} sort -i ${set}/hum/${feat}.hg38.split.bed > ${set}/hum/${feat}.hg38.sort.bed

# re-map: hum to dog
CrossMap.py bed ${hum2dog} ${set}/hum/${feat}.hg38.sort.bed ${set}/dog/${feat}.hg38-to-cf3.bed

# sort
${bedtools} sort -i ${set}/dog/${feat}.hg38-to-cf3.bed > ${set}/dog/${feat}.hg38-to-cf3.sort.bed

# merge
${bedtools} merge -i ${set}/dog/${feat}.hg38-to-cf3.sort.bed -c 10,11 -o distinct -d 10 > ${set}/dog/${feat}.hg38-to-cf3.sort.merge.bed

# back-map: dog to hum
CrossMap.py bed ${dog2hum} ${set}/dog/${feat}.hg38-to-cf3.sort.merge.bed ${set}/hum/${feat}.hg38-to-cf3-to-hg38.bed

# check reciprocal: intersect
${bedtools} sort -i ${set}/hum/${feat}.hg38-to-cf3-to-hg38.bed | ${bedtools} intersect -wa -a ${set}/hum/${feat}.hg38.sort.bed -b - | awk -F "\t" '{print $5}' | sort | uniq > ${set}/hum/${feat}.hg38-to-cf3-to-hg38.recip.txt

# keep reciprocal
grep -wf ${set}/hum/${feat}.hg38-to-cf3-to-hg38.recip.txt ${set}/dog/${feat}.hg38-to-cf3.sort.merge.bed > ${set}/dog/${feat}.hg38-to-cf3.sort.merge.recip.bed

# add flank and merge
${bedtools} sort -i ${set}/dog/${feat}.hg38-to-cf3.sort.merge.recip.bed | ${bedtools} slop -g cf3.chrom.sizes -b 1000 -i - | ${bedtools} merge -i - -c 4,5,6 -o distinct -d 10 > ${set}/dog/${feat}.hg38-to-cf3.sort.merge.recip.final.bed

# get matching
awk -F "\t" '{print $5}' ${set}/dog/${feat}.hg38-to-cf3.shared.bed > ${set}/dog/${feat}.hg38-to-cf3.shared.txt


# try with liftOver instead:
# map hum -> dog
/seq/vgb/dd/gwas/annotate/liftOver/liftOver \
-minMatch=0.5 \
/seq/vgb/dd/gwas/annotate/regions/${set}/hum/${feat}.hg38.genes-only.formatted.bed \
/seq/vgb/dd/gwas/annotate/liftOver/chainFiles/hg38_to_cf3.chain \
/seq/vgb/dd/gwas/annotate/regions/${set}/dog/${feat}.hg38-to-cf3.lO-mapped.bed \
/seq/vgb/dd/gwas/annotate/regions/${set}/dog/${feat}.hg38-to-cf3.lO-unmapped.bed

# map back
/seq/vgb/dd/gwas/annotate/liftOver/liftOver \
-minMatch=0.5 \
/seq/vgb/dd/gwas/annotate/regions/${set}/dog/${feat}.hg38-to-cf3.lO-mapped.bed \
/seq/vgb/dd/gwas/annotate/liftOver/chainFiles/cf3_to_hg38.chain \
/seq/vgb/dd/gwas/annotate/regions/${set}/hum/${feat}.hg38-to-cf3-to-hg38.lO-mapped.bed \
/seq/vgb/dd/gwas/annotate/regions/${set}/hum/${feat}.hg38-to-cf3-to-hg38.lO-unmapped.bed

# check reciprocal
/seq/vgb/software/bedtools/bedtools sort \
-i /seq/vgb/dd/gwas/annotate/regions/${set}/hum/${feat}.hg38-to-cf3-to-hg38.lO-mapped.bed | \
/seq/vgb/software/bedtools/bedtools intersect \
-wb -a /seq/vgb/dd/gwas/annotate/regions/${set}/hum/${feat}.hg38.genes-only.bed -b - > /seq/vgb/dd/gwas/annotate/regions/${set}/hum/${feat}.hg38-to-cf3-to-hg38.lO-mapped.recip.B.bed

/seq/vgb/software/bedtools/bedtools sort \
-i /seq/vgb/dd/gwas/annotate/regions/${set}/hum/${feat}.hg38-to-cf3-to-hg38.lO-mapped.bed | \
/seq/vgb/software/bedtools/bedtools intersect \
-wb -a - -b /seq/vgb/dd/gwas/annotate/regions/${set}/hum/${feat}.hg38.genes-only.bed > /seq/vgb/dd/gwas/annotate/regions/${set}/hum/${feat}.hg38-to-cf3-to-hg38.lO-mapped.recip.A.bed

# sort + merge canine coordinates


# map human coordinates to canine coordinates
CrossMap.py bed /seq/vgb/dd/gwas/annotate/liftOver/chainFiles/hg38_to_cf3.chain hum/${feat}.hg38.parsed.bed dog/${feat}.hg38-to-cf3.bed

# sort + merge canine coordinates
/seq/vgb/software/bedtools/bedtools sort -i dog/${feat}.hg38-to-cf3.bed | /seq/vgb/software/bedtools/bedtools merge -i - -c 4 -o distinct -d 10 > dog/${feat}.hg38-to-cf3.sort.merge.bed

# back-map canine coordinates to human coordinates
CrossMap.py bed /seq/vgb/dd/gwas/annotate/liftOver/chainFiles/cf3_to_hg38.chain dog/${feat}.hg38-to-cf3.sort.merge.bed hum/${feat}.hg38-to-cf3-to-hg38.bed

# get valid reciprocal human feature records
/seq/vgb/software/bedtools/bedtools sort -i hum/${feat}.hg38-to-cf3-to-hg38.bed | /seq/vgb/software/bedtools/bedtools intersect -wa -a hum/gencode.genes.hg38.genes-only.bed -b - > hum/${feat}.hg38-to-cf3-to-hg38.recip.bed

library(tidyverse)
dog=read_tsv("dog/${feat}.hg38-to-cf3.sort.merge.bed", col_names = F)
hum=read_tsv("hum/${feat}.hg38-to-cf3-to-hg38.recip.bed", col_names = F)
shared = merge((dog %>% unique() %>% group_by(X4) %>% summarize(chr=unique(X1),start=min(X2),end=max(X3)) %>% ungroup() %>% unique() %>% mutate(dog_len=end-start) %>% select(dog_chr=chr,dog_start=start,dog_end=end,dog_len,id=X4)),(hum %>% unique() %>% group_by(X4) %>% summarize(chr=unique(X1),start=min(X2),end=max(X3)) %>% ungroup() %>% unique() %>% mutate(hum_len=end-start) %>% select(hum_chr=chr,hum_start=start,hum_end=end,hum_len,id=X4)), by = "id")

library(tidyverse)
dog=read_tsv("dog/gencode.genes.hg38-to-cf3.sort.merge.bed",col_names=F)
dog = dog %>%
unique() %>%
group_by(X4) %>%
summarize(chr=unique(X1),start=min(X2),end=max(X3)) %>%
ungroup() %>%
unique() %>%
mutate(dog_len=end-start) %>%
select(dog_chr=chr,dog_start=start,dog_end=end,dog_len,id=X4) %>%
filter(!id %like% ",")

hum=read_tsv("hum/gencode.genes.hg38-to-cf3-to-hg38.recip.uniq.bed",col_names=F) %>% unique() %>% group_by(X4) %>% summarize(chr=unique(X1),start=min(X2),end=max(X3)) %>% ungroup() %>% unique() %>% mutate(hum_len=end-start) %>% select(hum_chr=chr,hum_start=start,hum_end=end,hum_len,id=X4)

/seq/vgb/software/bedtools/bedtools slop -g /seq/vgb/dd/gwas/annotate/regions/cf3.chrom.sizes -b 1000 -i - |




/seq/vgb/software/gcta/current --bfile geno/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465 --mBAT-combo <(tail -n+2 assoc/2022-11-20/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_phe-bq.146.mean-binary_dcov-datatype_qcov-age.hgt.bq.146_rel-cutoff-0.75.loco.mlma | awk -F "\t" 'OFS=FS {print $2,$4,$5,$6,$7,$8,$9,"3465"}') --mBAT-gene-list /seq/vgb/dd/gwas/annotate/regions/GENCODE/GENCODE.shared.hg38-to-cf3.reciprocal.dog.bed --out test
