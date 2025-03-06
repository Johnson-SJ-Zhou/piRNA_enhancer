#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate zbh_env

srna_dir="srna_antisense"
data_dir="/work/403spno1/thesis"
adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"


#trimming
cd $data_dir/dataset/$srna_dir/0_fastq/
for f in SRR*; do 
    cutadapt -a $adapter -q 20 -o $data_dir/dataset/$srna_dir/1_cutadapt/${f%.fastq}_trimmed.fastq $f
done


#alignment
cd $data_dir/dataset/$srna_dir/1_cutadapt/
for f in *; do 
    tophat2 -o $data_dir/dataset/$srna_dir/2_tophat/${file%.fastq} -p 16 -g 1 -N 0 --read-gap-length 0 --read-edit-dist 0 --read-realign-edit-dist 0 -G $data_dir/databases/GRCm39.gtf \
    --max-insertion-length 0 --max-deletion-length 0 --transcriptome-index=$db_dir/tophat_transcriptome_index/GRCm39 \
    $data_dir/databases/bt2l_index/GRCm39 $f
done


#bam to fastq
cd $data_dir/dataset/$srna_dir/2_tophat/
for file in *bam; do 
	bamToFastq -i $file -fq ${file%.bam}.fastq 
done


#Remove rRNA and tRNA
for file in *.fastq; do
	bowtie --norc -S -v 0 -p 16 --al $data_dir/dataset/$srna_dir/3_bowtie/${file%.fastq}_rrna_trna.fastq --un $data_dir/dataset/$srna_dir/3_bowtie/${file%.fastq}_rrna_trna_removed.fastq \
	$data_dir/databases/Rfam/Rfam_trna_rrna $file $data_dir/dataset/$srna_dir/3_bowtie/${file%.fastq}_rrna_trna.sam > $data_dir/dataset/$srna_dir/3_bowtie/${file%.fastq}_rRNAtRNA.log 2> $data_dir/dataset/$srna_dir/3_bowtie/${file%.fastq}_rRNAtRNA.err
done


#miRNA prediction and calculation
cd $data_dir/dataset/$srna_dir/4_mirdeep/
for file in $data_dir/dataset/$srna_dir/3_bowtie/*removed.fastq; do
	f1=`basename $file`
	mir_dir=$data_dir/databases/mirbase
	geno_dir=$data_dir/databases/ensembl/mirdeep2_genome
	mapper.pl $file -e -h -j -l 18 -m -p $geno_dir/mm39_genome_mirdeep -s ${f1%.fastq}_reads_collapsed.fa -t ${f1%.fastq}_reads_collapsed_vs_genome.arf -v -o 16
	miRDeep2.pl ${f1%.fastq}_reads_collapsed.fa  $geno_dir/mm39_genome_mirdeep.fa ${f1%.fastq}_reads_collapsed_vs_genome.arf \
	$mir_dir/mature_mmu.fa $mir_dir/mature_hsa_rno.fa  $mir_dir/hairpin_mmu.fa  -t mmu 2> ${f1%.fastq}_report.log
	quantifier.pl -p $mir_dir/hairpin_mmu.fa -m $mir_dir/mature_mmu.fa -r ${f1%.fastq}_reads_collapsed.fa
done


#miRNA precursor prediction
cd mirdeep_runs
for file in * ; do
	cat $data_dir/dataset/$srna_dir/4_mirdeep/mirdeep_runs/$file/output.mrd | egrep "^pri_seq" | cut -c31- | tr "u" "t" | tr "atcg" "ATCG" | perl $data_dir/program/seq2fasta.pl > $data_dir/dataset/$srna_dir/4_mirdeep/mirdeep_runs/$file/$file"_miRDeepPredictedPrecursors.fa"
	name=`cat $data_dir/dataset/$srna_dir/4_mirdeep/mirdeep_runs/$file/$file"_parameters" | egrep "^file_reads_vs_genome" | cut -c22- `
	cat $data_dir/dataset/$srna_dir/4_mirdeep/mirdeep_runs/$file/$file"_miRDeepPredictedPrecursors.fa" > $data_dir/dataset/$srna_dir/4_mirdeep/${name%_reads_collapsed_vs_genome.arf}"_miRDeepPredictedPrecursors.fa" 
done


##filter_miRNA
cd $data_dir/dataset/$srna_dir/3_bowtie/
for file in *_rrna_trna_removed.fastq; do
	bowtie -norc -S -v 0 -p 16 --al ${file%_rrna_trna_removed.fastq}_miRBase.fastq --un ${file%_rrna_trna_removed.fastq}_miRBase_removed.fastq ${db_dir%mm39/Gencode}mirbase/miRBase_pre \
	$file ${file%_rrna_trna_removed.fastq}_miRBase.sam > ${file%_rrna_trna_removed.fastq}_miRBase_removed.log 2> ${file%_rrna_trna_removed.fastq}_miRBase_removed.err
done

cd $data_dir/dataset/$srna_dir/4_mirdeep/
mkdir bt_index
for file in *_miRDeepPredictedPrecursors.fa; do
	f1=${file%_rrna_trna_removed_miRDeepPredictedPrecursors.fa}
	bowtie-build $file $data_dir/dataset/$srna_dir/4_mirdeep/bt_index/${file%.fa}
	bowtie -norc -S -v 0 -p 16 --al $data_dir/dataset/$srna_dir/3_bowtie/$f1"_miRDeep2.fastq" --un $data_dir/dataset/$srna_dir/3_bowtie/$f1"_miRDeep2_removed.fastq" $data_dir/dataset/$srna_dir/4_mirdeep/bt_index/${file%.fa} $data_dir/dataset/$srna_dir/3_bowtie/$f1"_miRBase_removed.fastq" \
	$data_dir/dataset/$srna_dir/3_bowtie/$f1"_miRDeep2.sam" > $data_dir/dataset/$srna_dir/3_bowtie/$f1"_miRDeep2_removed.log" 2> $data_dir/dataset/$srna_dir/3_bowtie/$f1"_miRDeep2_removed.err"
done


#size select
cd $data_dir/dataset/$srna_dir/3_bowtie/
for file in *miRDeep2_removed.fastq; do
	cat $file | perl $data_dir/program/filterFastqByLen.pl 24 34 > $data_dir/dataset/$srna_dir/5_size_select/${file%_miRDeep2_removed.fastq}_filteredByLen.fastq 2> $data_dir/dataset/$srna_dir/5_len_filter/${file%_miRDeep2_removed.fastq}_filteredByLen.err
done


##align_to_genome
cd $data_dir/dataset/$srna_dir/5_size_select/
for file in *fastq; do
	tophat2 -o $data_dir/dataset/$srna_dir/2_tophat/${file%.fastq} -p 16 -g 1 -N 0 \
	--read-gap-length 0 --read-edit-dist 0 --read-realign-edit-dist 0 -G $data_dir/databases/GRCm39.gtf \
	--max-insertion-length 0 --max-deletion-length 0 --transcriptome-index=$db_dir/tophat_transcriptome_index/GRCm39 \
	$data_dir/databases/bt2l_index/GRCm39 $file
done


##bam to bed
cd $data_dir/dataset/$srna_dir/2_tophat/
for f in *_filteredByLen.bam; do
  bedtools bamtobed -i $data_dir/dataset/$srna_dir/2_tophat/$f/accepted_hits.bam > $data_dir/dataset/$srna_dir/3_bed/$f".bed"
done


#intersect SE
cd $data_dir/dataset/$srna_dir/5_size_select/
for f in *_filteredByLen.bed; do
  cd $data_dir/databases/SEdb/mm39_SE/merged
  for i in *.bed; do
    bedtools intersect -u -f 0.8  -a $data_dir/dataset/$srna_dir/5_size_select/$f -b $data_dir/databases/SEdb/mm39_SE/merged/$i > $data_dir/dataset/$srna_dir/6_intersect_SE/${f%.bed}"_"${i%.bed}".bed"
  done
done

wc -l $data_dir/dataset/$srna_dir/6_intersect_SE/* > ~/srna_seq_antisense"_6_intersect_SE_merged.txt"


#intersect TE
cd $data_dir/dataset/$srna_dir/5_size_select/
for f in *_filteredByLen.bed; do
  cd $data_dir/databases/SEdb/mm39_TE/merged
  for i in *.bed; do
    bedtools intersect -u -f 0.8  -a $data_dir/dataset/$srna_dir/5_size_select/$f -b $data_dir/databases/SEdb/mm39_TE/merged/$i > $data_dir/dataset/$srna_dir/6_intersect_TE/${f%.bed}"_"${i%.bed}".bed" 
  done
done

wc -l $data_dir/dataset/$srna_dir/6_intersect_TE/* > ~/srna_seq_antisense"_6_intersect_TE_merged.txt"


#intersect meiosis mitosis SE
cd $data_dir/dataset/$srna_dir/5_size_select/
for f in *_filteredByLen.bed; do
  cd $data_dir/databases/meiosis_mitosis
  for i in *.bed; do
    bedtools intersect -u -f 0.8  -a $data_dir/dataset/$srna_dir/5_size_select/$f -b $data_dir/databases/meiosis_mitosis/$i > $data_dir/dataset/$srna_dir/6_intersect_meiosis_mitosis_SE/${f%.bed}"_"${i%.bed}".bed"
  done
done

wc -l $data_dir/dataset/$srna_dir/6_intersect_meiosis_mitosis_SE/* > ~/srna_seq_antisense"_6_intersect_mm_SE_merged.txt"


#enhancer pirna abundance count
cd $data_dir/dataset/$srna_dir/5_size_select/
for f in *_filteredByLen.bed; do
  bedtools intersect -F 0.8 -a $data_dir/databases/SEdb/mm39_SE/count/all.sorted.bed -b $data_dir/dataset/$srna_dir/5_size_select/$f -c -bed > $data_dir/dataset/$srna_dir/6_enhancer_counts_SE/$f
  bedtools intersect -F 0.8 -a $data_dir/databases/SEdb/mm39_TE/count/all.sorted.bed -b $data_dir/dataset/$srna_dir/5_size_select/$f -c -bed > $data_dir/dataset/$srna_dir/6_enhancer_counts_TE/$f
done


#intersect erna
cd $data_dir/dataset/$srna_dir/5_size_select/
for f in *_filteredByLen.bed; do
  cd $data_dir/databases/eRNAbase/mm39/merged
  for i in *.bed; do
    bedtools intersect -u -f 0.8 -s -a $data_dir/dataset/$srna_dir/5_size_select/$f -b $data_dir/databases/eRNAbase/mm39/merged/$i > $data_dir/dataset/$srna_dir/7_intersect_erna/${f%.bed}"_"${i%.bed}".bed"
  done
done
wc -l $data_dir/dataset/$srna_dir/7_intersect_erna/* > ~/srna_seq_antisense"_7_intersect_erna.txt"


#enhancer pirna abundance count
cd $data_dir/dataset/$srna_dir/5_size_select/
for f in *_filteredByLen.bed; do
  bedtools intersect -F 0.8 -s -a $data_dir/databases/eRNAbase/mm39/all_peaks_gene.bed -b $data_dir/dataset/$srna_dir/5_size_select/$f -c -bed > $data_dir/dataset/$srna_dir/7_enhancer_counts_erna/$f
done


