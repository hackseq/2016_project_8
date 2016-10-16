from __future__ import print_function # Python 2.x
# import pysam
import vcf
import tenkit
import pyfasta
import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam
import subprocess


def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    stdout_lines = iter(popen.stdout.readline, "")
    for stdout_line in stdout_lines:
        yield stdout_line

    popen.stdout.close()
    return_code = popen.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, cmd)


def run_freebayes(reference_path, out_vcf, bed_path, bam_path):
    """
    Use freebayes to create a VCF file
    :param reference_path: path to the reference file
    :param out_vcf: path to the output vcf file
    :param bed_path: path to the bed file
    :param bam_path: path to the bam file
    :return: VCF file as output on disk
    """
    command = ['freebayes', '-f', reference_path,
               '-0', '-C', '3', '-F', '0.03', '--pooled-continuous', '--pooled-discrete',
               '--min-coverage', '10', '-v', out_vcf, '-t', bed_path, '-b', bam_path]
    for output in execute(command):
        print(output)

# Open a BAM file
# pct0 = "/hackseq/HCC1954_Exome_Data_for_HackSeq/HCC1954_0pct_phased_possorted.bam"
# p0_bam = pysam.Samfile(pct0)


# Open a VCF file
# pct0v = "/hackseq/HCC1954_Exome_Data_for_HackSeq/HCC1954_0pct_snpindel.vcf.gz"
# p0_vcf = vcf.Reader(open(pct0v))


# Open the reference
# fa = pyfasta.Fasta("/hackseq/hg19/refdata-hg19-2.1.0/fasta/genome.fa")

# Get the table h1/h2 ref/alt counts for a variant in vcf_record in the given bam file.
def get_counts_for_record(vcf_record, bam):
    alleles = tk_io.get_record_alt_alleles(rec)
    ref = tk_io.get_record_ref(rec)
    r = get_allele_read_info(rec.CHROM, rec.POS, ref, alleles, 30, p0_bam, fa)

    
# The counts of each allele on each haplotype
def get_allele_read_info(chrom, pos, ref, alt_alleles, min_mapq, bam, 
                         reference_pyfasta, max_reads=2000, match = 1, 
                         mismatch = -4, gap_open = -6, gap_extend = -1):
    all_alleles = [ref] + alt_alleles
    
    counts = {'h1': [0 for x in all_alleles], 
              'h2': [0 for x in all_alleles], 
              'un': [0 for x in all_alleles]}
    
    num_reads = 0
    qnames = set()
    
    # Fetch all the reads
    for read in bam.fetch(chrom, pos, pos + 1):
        num_reads += 1
        
        # Only let each read count once, and filter duplicated reads
        if read.qname in qnames:
            continue
        qnames.add(read.qname)
        if read.is_duplicate:
            continue
        if num_reads > max_reads:
            break

        is_indel_variant = False
        for allele in alt_alleles:
            if len(allele) != len(ref):
                is_indel_variant = True
    
        # Get haplotype assigned to read
        hap = tk_io.get_read_haplotype(read)
        if hap == 1:
            hap = "h1"
        elif hap == 2:
            hap = "h2"
        elif hap == None:
            hap = "un"
        else:
            print("unknown hap: %s" % str(hap))
                
        # This aligns the read sequence to both alleles to avoid alignment artifacts
        allele_index_in_read = tk_bam.read_contains_allele_sw(ref, all_alleles, pos, 
                                                              read, reference_pyfasta[chrom], 
                                                              match = match, mismatch = mismatch, 
                                                              gap_open = gap_open, gap_extend = gap_extend)
        for (allele_index, allele) in enumerate(all_alleles):
            if allele_index == allele_index_in_read:
                    
                if read.mapq >= min_mapq:
                    counts[hap][allele_index] += 1
                    
    counts_reformat = { k: {'ref': v[0], 'alt': v[1]} for (k,v) in counts.items() }               
    return counts_reformat

# rec = p0_vcf.next()
# alleles = tk_io.get_record_alt_alleles(rec)
# ref = tk_io.get_record_ref(rec)
# r = get_allele_read_info(rec.CHROM, rec.POS, ref, alleles, 30, p0_bam, fa)
# r = get_allele_read_info('chr17', 41258433, 'G', ['T'], 30, p0_bam, fa)
# get_ipython().magic(u'pinfo2 tk_io.get_read_haplotype')


if __name__ == '__main__':
    run_freebayes('/hackseq/hg19/refdata-hg19-2.1.0/fasta/genome.fa', '/hackseq/team8somatic/super-vcf-out2.vcf',
                  '/hackseq/team8somatic/bed_files/exome_chr22.bed',
                  '/hackseq/HCC1954_Exome_Data_for_HackSeq/HCC1954_0pct_phased_possorted.bam')