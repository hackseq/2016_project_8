# Links!!!!


Datasets:
[http://52.89.192.184/hackseq/HCC1954_Exome_Data_for_HackSeq/](http://52.89.192.184/hackseq/HCC1954_Exome_Data_for_HackSeq/)

IGV:
[http://software.broadinstitute.org/software/igv/download_snapshot](http://software.broadinstitute.org/software/igv/download_snapshot)

IPython Notebook:
[http://ec2-52-89-192-184.us-west-2.compute.amazonaws.com/](http://ec2-52-89-192-184.us-west-2.compute.amazonaws.com/)


BAM documentation:
http://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam

VCF documentation:
http://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/vcf


Getting the longranger 'environment'. On EC2 instance:
source /hackseq/longranger-2.1.1/sourceme.bash


# Pipeline Planning

## run_freebayes(ref, bam, bed, out_vcf) - Pat

Runs freebayes to shortlist positions specified by a bed file.
Processes the resulting VCF file into a simpler format with chrom, pos, ref, alt information.

## count_position(pos, bam, ref) - Aman, Peng

Taking one position and counts the haplotype 1 and haplotype 2 read counts.

## model_positions(position_counts)

Runs the statistical model on each position's read counts
For each position, determines the max likelihood model.

## call_all(position_list, bam, ref) - Aman

Wrapper around count_position() which loops it on all the the shortlisted positions.

## write_output

Outputs into tabular format for downstream statistical analyses and sanity checks.
