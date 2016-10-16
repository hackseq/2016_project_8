
# coding: utf-8

# In[48]:

# tour of pysam
import pysam
import vcf
import tenkit
import pyfasta
import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam


# In[66]:

# Open a BAM file
pct0 = "/hackseq/HCC1954_Exome_Data_for_HackSeq/HCC1954_0pct_phased_possorted.bam"
p0_bam = pysam.Samfile(pct0)


# In[67]:

# Open a VCF file
pct0v = "/hackseq/HCC1954_Exome_Data_for_HackSeq/HCC1954_0pct_snpindel.vcf.gz"
p0_vcf = vcf.Reader(open(pct0v))


# In[51]:

# Open the reference
fa = pyfasta.Fasta("/hackseq/hg19/refdata-hg19-2.1.0/fasta/genome.fa")


# In[1]:

#
# Get the table h1/h2 ref/alt counts for a variant in vcf_record in the given bam file. 
#
def get_counts_for_record(vcf_record, bam):
    alleles = tk_io.get_record_alt_alleles(rec)
    ref = tk_io.get_record_ref(rec)
    r = get_allele_read_info(rec.CHROM, rec.POS, ref, alleles, 30, p0_bam, fa)
    r
    
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
            print "unknown hap: %s" % str(hap)
                
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


# In[59]:

rec = p0_vcf.next()


# In[60]:

rec.POS


# In[61]:

alleles = tk_io.get_record_alt_alleles(rec)
ref = tk_io.get_record_ref(rec)
r = get_allele_read_info(rec.CHROM, rec.POS, ref, alleles, 30, p0_bam, fa)
r


# In[68]:

r = get_allele_read_info('chr17', 41258433, 'G', ['T'], 30, p0_bam, fa)
r


# In[24]:

get_ipython().magic(u'pinfo2 tk_io.get_read_haplotype')


# In[ ]:



