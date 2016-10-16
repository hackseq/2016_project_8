"""
Hopefully a main file to run everything.
"""

from count import run_freebayes


if __name__ == '__main__':
    run_freebayes('/hackseq/hg19/refdata-hg19-2.1.0/fasta/genome.fa',
                  '/hackseq/team8somatic/super-vcf-out2.vcf',
                  '/hackseq/team8somatic/bed_files/exome_chr22.bed',
                  '/hackseq/HCC1954_Exome_Data_for_HackSeq/HCC1954_0pct_phased_possorted.bam')