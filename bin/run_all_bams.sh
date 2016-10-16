vcf = $1

python main.py /hackseq/hg19/refdata-hg19-2.1.0/fasta/genome.fa /hackseq/HCC1954_Exome_Data_for_HackSeq/749709.vcf.gz /hackseq/HCC1954_Exome_Data_for_HackSeq/HCC1954_0pct_phased_possorted.bam /hackseq/team8somatic/count_table_0pct.csv &

python main.py /hackseq/hg19/refdata-hg19-2.1.0/fasta/genome.fa /hackseq/HCC1954_Exome_Data_for_HackSeq/749709.vcf.gz /hackseq/HCC1954_Exome_Data_for_HackSeq/HCC1954_25pct_phased_possorted.bam /hackseq/team8somatic/count_table_25pct.csv &

python main.py /hackseq/hg19/refdata-hg19-2.1.0/fasta/genome.fa /hackseq/HCC1954_Exome_Data_for_HackSeq/749709.vcf.gz /hackseq/HCC1954_Exome_Data_for_HackSeq/HCC1954_50pct_phased_possorted.bam /hackseq/team8somatic/count_table_50pct.csv &

python main.py /hackseq/hg19/refdata-hg19-2.1.0/fasta/genome.fa /hackseq/HCC1954_Exome_Data_for_HackSeq/749709.vcf.gz /hackseq/HCC1954_Exome_Data_for_HackSeq/HCC1954_75pct_phased_possorted.bam /hackseq/team8somatic/count_table_75pct.csv &

python main.py /hackseq/hg19/refdata-hg19-2.1.0/fasta/genome.fa /hackseq/HCC1954_Exome_Data_for_HackSeq/749709.vcf.gz /hackseq/HCC1954_Exome_Data_for_HackSeq/HCC1954_100pct_phased_possorted.bam /hackseq/team8somatic/count_table_100pct.csv &
