#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate zbh_env

IP_dir="piwi_IP_pirbase_antisense"
adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
data_dir="/work/403spno1/thesis"


#trimming
cd $data_dir/dataset/chip/0_fastq/
for f in SRR*; do 
    cutadapt -a $adapter -q 20 -o $data_dir/dataset/chip/1_cutadapt/${f%.fastq}_trimmed.fastq $f
done


#size select
cd $data_dir/dataset/$IP_dir/1_cutadapt/
for f in SRR*; do 
    cat $f | perl $data_dir/scripts/filterFastqByLen.pl 24 34 > $data_dir/dataset/$IP_dir/2_size_select/${f%_trimmed.fastq}_filteredByLen.fastq 2> $data_dir/dataset/$IP_dir/2_size_select/${f%_trimmed.fastq}_filteredByLen.err
done

#alignment
cd $data_dir/dataset/$IP_dir/2_size_select/
for f in *; do 
    tophat2 -o $data_dir/2_tophat/${f%.fastq} -p 16 -g 1 -N 0 --read-gap-length 0 --read-edit-dist 0 --read-realign-edit-dist 0 -G $data_dir/databases/GRCm39.gtf \
    --max-insertion-length 0 --max-deletion-length 0 --transcriptome-index=$db_dir/tophat_transcriptome_index/GRCm39 \
    $data_dir/databases/bt2l_index/GRCm39 $f
done


#bam to bed
cd $data_dir/dataset/$IP_dir/2_tophat/
for f in *; do
  bedtools bamtobed -i $data_dir/dataset/$IP_dir/2_tophat/$f/accepted_hits.bam > $data_dir/dataset/$IP_dir/3_bed/$f".bed"
done 


#convert chrom names
cd $data_dir/dataset/$IP_dir/3_bed/
for f in *.bed; do
  python $data_dir/scripts/07_convert_chrom.py $f ${f%.bed}_temp.bed
  awk '{ s=$3-$2 } s >= 24 && s <= 34 { print }' ${f%.bed}_temp.bed > ${f%.bed}_corrected.bed
  rm ${f%.bed}_temp.bed
done

wc -l $data_dir/dataset/$IP_dir/3_bed/*corrected.bed > ~/$IP_dir"_3_bed_total.txt"


#intersect SE
cd $data_dir/dataset/$IP_dir/3_bed/
for f in *no_mrna_lncrna.bed; do
  cd $data_dir/databases/SEdb/mm39_SE/merged
  for i in *.bed; do
    bedtools intersect -u -f 0.8 -s  -a $data_dir/dataset/$IP_dir/3_bed/$f -b $data_dir/databases/SEdb/mm39_SE/merged/$i > $data_dir/dataset/$IP_dir/5_intersect_SE/${f%.bed}"_"${i%.bed}".bed"
    bedtools intersect -u -f 0.8  -a $data_dir/dataset/$IP_dir/3_bed/$f -b $data_dir/databases/SEdb/mm39_SE/merged_unique/$i > $data_dir/dataset/$IP_dir/5_intersect_SE/merged_unique/${f%.bed}"_"${i%.bed}".bed"
  done
done

wc -l $data_dir/dataset/$IP_dir/5_intersect_SE/merged/* > ~/$IP_dir"_5_intersect_SE_merged.txt"
wc -l $data_dir/dataset/$IP_dir/5_intersect_SE/merged_unique/* > ~/$IP_dir"_5_intersect_SE_merged_unique.txt"


#intersect TE
cd $data_dir/dataset/$IP_dir/3_bed/
for f in *no_mrna_lncrna.bed; do
  cd $data_dir/databases/SEdb/mm39_TE/merged
  for i in *.bed; do
    bedtools intersect -u -f 0.8 -s -a $data_dir/dataset/$IP_dir/3_bed/$f -b $data_dir/databases/SEdb/mm39_TE/merged/$i > $data_dir/dataset/$IP_dir/6_intersect_TE/${f%.bed}"_"${i%.bed}".bed" 
    bedtools intersect -u -f 0.8  -a $data_dir/dataset/$IP_dir/3_bed/$f -b $data_dir/databases/SEdb/mm39_TE/merged_unique/$i > $data_dir/dataset/$IP_dir/6_intersect_TE/merged_unique/${f%.bed}"_"${i%.bed}".bed"
  done
done

wc -l $data_dir/dataset/$IP_dir/6_intersect_TE/merged/* > ~/$IP_dir"_6_intersect_TE_merged.txt"
wc -l $data_dir/dataset/$IP_dir/6_intersect_TE/merged_unique/* > ~/$IP_dir"_6_intersect_TE_merged_unique.txt"


#intersect meiosis mitosis SE
cd $data_dir/dataset/$IP_dir/3_bed/
for f in *corrected.bed; do
  cd $data_dir/databases/meiosis_mitosis
  for i in *.bed; do
    bedtools intersect -u -f 0.8  -a $data_dir/dataset/$IP_dir/3_bed/$f -b $data_dir/databases/meiosis_mitosis/$i > $data_dir/dataset/$IP_dir/5_intersect_meiosis_mitosis_SE/${f%.bed}"_"${i%.bed}".bed"
  done
done

wc -l $data_dir/dataset/$IP_dir/5_intersect_meiosis_mitosis_SE/* > ~/$IP_dir"_5_intersect_mm_SE_merged.txt"


#enhancer pirna abundance count
cd $data_dir/dataset/$IP_dir/3_bed/
for f in *corrected.bed; do
  bedtools intersect -F 0.8 -a $data_dir/databases/SEdb/mm39_SE/count/all.sorted.bed -b $data_dir/dataset/$IP_dir/3_bed/$f -c -bed > $data_dir/dataset/$IP_dir/5_enhancer_counts_SE/$f
  bedtools intersect -F 0.8 -a $data_dir/databases/SEdb/mm39_TE/count/all.sorted.bed -b $data_dir/dataset/$IP_dir/3_bed/$f -c -bed > $data_dir/dataset/$IP_dir/6_enhancer_counts_TE/$f
done

cd $data_dir/dataset/$IP_dir/3_bed/
for f in *corrected.bed; do bedtools intersect -F 0.8 -a $data_dir/databases/meiosis_mitosis/meiotic_specific_SE_mm39.bed -b $data_dir/dataset/$IP_dir/3_bed/$f -c -bed > $data_dir/dataset/$IP_dir/5_enhancer_counts_MM/$f; done


#intersect erna
cd $data_dir/dataset/$IP_dir/3_bed/
for f in *corrected.bed; do
  cd $data_dir/databases/eRNAbase/mm39/merged
  for i in *corrected.bed; do
    bedtools intersect -u -f 0.8 -s -a $data_dir/dataset/$IP_dir/3_bed/$f -b $data_dir/databases/eRNAbase/mm39/merged/$i > $data_dir/dataset/$IP_dir/7_intersect_erna/${f%.bed}"_"${i%.bed}".bed"
  done
done

wc -l $data_dir/dataset/$IP_dir/7_intersect_erna/* > ~/$IP_dir"_7_intersect_erna.txt"

#erna pirna abundance count
cd $data_dir/dataset/$IP_dir/3_bed/
for f in *corrected.bed; do
  bedtools intersect -F -0.8 -s -a $data_dir/databases/eRNAbase/mm39/all_peaks_gene.bed -b $data_dir/dataset/$IP_dir/3_bed/$f -c -bed > $data_dir/dataset/$IP_dir/7_enhancer_counts_erna/$f
done

