#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate zbh_env

srna_dir="srna_antisense_rnacentral"
data_dir="$data_dir"
adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"


#trimming
cd $data_dir/0_fastq/
for f in SRR*; do 
    cutadapt -a $adapter -q 20 -o $data_dir/1_cutadapt/${f%.fastq}_trimmed.fastq $f
done


#alignment
cd $data_dir/1_cutadapt/
for f in *; do 
    tophat2 -o $data_dir/2_tophat/${file%.fastq} -p 16 -g 1 -N 0 --read-gap-length 0 --read-edit-dist 0 --read-realign-edit-dist 0 -G $data_dir/databases/GRCm39.gtf \
    --max-insertion-length 0 --max-deletion-length 0 --transcriptome-index=$db_dir/tophat_transcriptome_index/GRCm39 \
    $data_dir/databases/bt2l_index/GRCm39 $f
done


#bam to bed
cd $data_dir/dataset/$srna_dir/2_tophat/
for f in *; do
  bedtools bamtobed -i $data_dir/dataset/$srna_dir/2_tophat/$f/accepted_hits.bam > $data_dir/dataset/$srna_dir/3_bed/$f".bed"
done


#remove rrna trna mirna
cd $data_dir/dataset/$srna_dir/3_bed/
for f in *.bed; do
  bedtools intersect -a $file -b $data_dir/databases/RNAcentral/rrna_trna_mirna.bed -f 0.8 -s -v > $data_dir/dataset/$srna_dir/4_remove_rna/${f%.bed}"_removed_rna.bed"
done


#size select
cd $data_dir/dataset/$srna_dir/4_remove_rna/
for f in *removed_rna.bed; do
  awk '{ s=$3-$2 } s >= 24 && s <= 34 { print }' $data_dir/dataset/$srna_dir/4_remove_rna/$f > $data_dir/dataset/$srna_dir/5_size_select/${f%_removed_rna.bed}_ss.bed
done

wc -l $data_dir/dataset/$srna_dir/5_size_select/*_ss.bed > ~/srna_seq_antisense"_total.txt"


#intersect SE
cd $data_dir/dataset/$srna_dir/5_size_select/
for f in *_ss.bed; do
  cd $data_dir/databases/SEdb/mm39_SE/merged
  for i in *.bed; do
    bedtools intersect -u -f 0.8 -a $data_dir/dataset/$srna_dir/5_size_select/$f -b $data_dir/databases/SEdb/mm39_SE/merged/$i > $data_dir/dataset/$srna_dir/5_intersect_SE/${f%.bed}"_"${i%.bed}".bed"
  done
done

wc -l $data_dir/dataset/$srna_dir/5_intersect_SE/*.bed > ~/srna_seq_antisense"_5_intersect_SE_merged.txt"


#intersect TE
cd $data_dir/dataset/$srna_dir/5_size_select/
for f in *ss.bed; do
  cd $data_dir/databases/SEdb/mm39_TE/merged
  for i in *.bed; do
    bedtools intersect -u -f 0.8 -a $data_dir/dataset/$srna_dir/5_size_select/$f -b $data_dir/databases/SEdb/mm39_TE/merged/$i > $data_dir/dataset/$srna_dir/6_intersect_TE/${f%.bed}"_"${i%.bed}".bed" 
  done
done

wc -l $data_dir/dataset/$srna_dir/6_intersect_TE/*.bed > ~/srna_seq_antisense"_6_intersect_TE_merged.txt"


#intersect meiosis mitosis SE
cd $data_dir/dataset/$srna_dir/5_size_select/
for f in *ss.bed; do
  cd $data_dir/databases/meiosis_mitosis
  for i in *.bed; do
    bedtools intersect -u -f 0.8  -a $data_dir/dataset/$srna_dir/5_size_select/$f -b $data_dir/databases/meiosis_mitosis/$i > $data_dir/dataset/$srna_dir/5_intersect_meiosis_mitosis_SE/${f%.bed}"_"${i%.bed}".bed"
  done
done

wc -l $data_dir/dataset/$srna_dir/5_intersect_meiosis_mitosis_SE/*.bed > ~/srna_seq_antisense"_5_intersect_mm_SE_merged.txt"


#enhancer pirna abundance count
cd $data_dir/dataset/$srna_dir/5_size_select/
for f in *ss.bed; do
  bedtools intersect -F 0.8 -a $data_dir/databases/SEdb/mm39_SE/count/all.sorted.bed -b $data_dir/dataset/$srna_dir/5_size_select/$f -c -bed > $data_dir/dataset/$srna_dir/5_enhancer_counts_SE/$f
  bedtools intersect -F 0.8 -a $data_dir/databases/SEdb/mm39_TE/count/all.sorted.bed -b $data_dir/dataset/$srna_dir/5_size_select/$f -c -bed > $data_dir/dataset/$srna_dir/6_enhancer_counts_TE/$f
done


#intersect erna
cd $data_dir/dataset/$srna_dir/5_size_select/
for f in *ss.bed; do
  cd $data_dir/databases/eRNAbase/mm39/merged
  for i in *.bed; do
    bedtools intersect -u -f 0.8 -s -a $data_dir/dataset/$srna_dir/5_size_select/$f -b $data_dir/databases/eRNAbase/mm39/merged/$i > $data_dir/dataset/$srna_dir/7_intersect_erna/${f%.bed}"_"${i%.bed}".bed"
  done
done
wc -l $data_dir/dataset/$srna_dir/7_intersect_erna/*.bed > ~/srna_seq_antisense"_7_intersect_erna.txt"


#enhancer pirna abundance count
cd $data_dir/dataset/$srna_dir/5_size_select/
for f in *ss.bed; do
  bedtools intersect -F 0.8 -s -a $data_dir/databases/eRNAbase/mm39/all_peaks_gene.bed -b $data_dir/dataset/$srna_dir/5_size_select/$f -c -bed > $data_dir/dataset/$srna_dir/7_enhancer_counts_erna/$f
done
