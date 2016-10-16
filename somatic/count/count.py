import pysam
import vcf
import tenkit
import pyfasta
import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam


def get_counts_for_record(vcf_rec, bam, fa):
    """
    Get the table h1/h2 ref/alt counts for a variant in vcf_record in the given bam file.
    :param vcf_rec:
    :param bam:
    :param fa:
    :return:
    """
    alleles = tk_io.get_record_alt_alleles(vcf_rec)
    ref = tk_io.get_record_ref(vcf_rec)
    print(vcf_rec.CHROM)
    r = get_allele_read_info(vcf_rec.CHROM, vcf_rec.POS, ref, alleles, 30, bam, fa)
    return r

def get_all_counts(px_vcf, bam, fa):
    for vcf in px_vcf:
        record = get_counts_for_record(vcf_rec, px_bam, fa)

    
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

    # import pdb; pdb.set_trace()

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
                                                              gap_open = gap_open,
                                                              gap_extend = gap_extend)
        for (allele_index, allele) in enumerate(all_alleles):
            if allele_index == allele_index_in_read:
                    
                if read.mapq >= min_mapq:
                    counts[hap][allele_index] += 1
                    
    counts_reformat = { k: {'ref': v[0], 'alt': v[1]} for (k,v) in counts.items() }               
    return counts_reformat


if __name__ == '__main__':
    # Open the reference
    ref = '/hackseq/hg19/refdata-hg19-2.1.0/fasta/genome.fa'
    fa = pyfasta.Fasta(ref)
    # Open a VCF file
    vcf_path = '/hackseq/team8somatic/super-vcf-out2.vcf'
    px_vcf = vcf.Reader(open(vcf_path))
    vcf_rec = px_vcf.next()
    # Open a BAM file
    bam = '/hackseq/HCC1954_Exome_Data_for_HackSeq/HCC1954_0pct_phased_possorted.bam'
    px_bam = pysam.Samfile(bam)
    # get_counts_for_record(vcf_rec, px_bam, fa)
    get_all_counts(px_vcf, px_bam, fa)