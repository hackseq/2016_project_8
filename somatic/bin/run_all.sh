#!/bin/bash
#

vcf=$1
name=$2
hcc_dir=/hackseq/HCC1954_Exome_Data_for_HackSeq
ref=/hackseq/hg19/refdata-hg19-2.1.0/fasta/genome.fa
count_table_dir=/hackseq/team8somatic/count_table
bed_option=--bed=/hackseq/team8somatic/20141020.strict_mask.whole_genome.bed

mkdir $count_table_dir/$name

for bam in 0 25 50 75 100
do
    python ../count/count.py $bed_option $ref $vcf $hcc_dir/HCC1954_${bam}pct_phased_possorted.bam $count_table_dir/$name/count_table_${bam}pct_$name.csv &
done

#python ../count/count.py $bed_option $ref $vcf $hcc_dir/HCC1954_0pct_phased_possorted.bam $count_table_dir/$name/count_table_0pct_$name.csv &
#
#python ../count/count.py $bed_option $ref $vcf $hcc_dir/HCC1954_25pct_phased_possorted.bam $count_table_dir/$name/count_table_25pct_$name.csv &
#
#python ../count/count.py $bed_option $ref $vcf $hcc_dir/HCC1954_50pct_phased_possorted.bam $count_table_dir/$name/count_table_50pct_$name.csv &
#
#python ../count/count.py $bed_option $ref $vcf $hcc_dir/HCC1954_75pct_phased_possorted.bam $count_table_dir/$name/count_table_75pct_$name.csv &
#
#python ../count/count.py $bed_option $ref $vcf $hcc_dir/HCC1954_100pct_phased_possorted.bam $count_table_dir/$name/count_table_100pct_$name.csv &