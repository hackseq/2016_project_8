vcf=$1
name=$2
hcc_dir=/hackseq/HCC1954_Exome_Data_for_HackSeq
ref=/hackseq/hg19/refdata-hg19-2.1.0/fasta/genome.fa


python ../count/count.py $ref $vcf $hcc_dir/HCC1954_0pct_phased_possorted.bam /hackseq/team8somatic/count_table/withrefalt/count_table_0pct_$name.csv &

python ../count/count.py $ref $vcf $hcc_dir/HCC1954_25pct_phased_possorted.bam /hackseq/team8somatic/count_table/withrefalt/count_table_25pct_$name.csv &

python ../count/count.py $ref $vcf $hcc_dir/HCC1954_50pct_phased_possorted.bam /hackseq/team8somatic/count_table/withrefalt/count_table_50pct_$name.csv &

python ../count/count.py $ref $vcf $hcc_dir/HCC1954_75pct_phased_possorted.bam /hackseq/team8somatic/count_table/withrefalt/count_table_75pct_$name.csv &

python ../count/count.py $ref $vcf $hcc_dir/HCC1954_100pct_phased_possorted.bam /hackseq/team8somatic/count_table/withrefalt/count_table_100pct_$name.csv &