#!/usr/bin/env python
#
#
import pysam
import vcf
import tenkit
import pyfasta
import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam
import pandas
import os
import sys
import docopt


def error(msg):
    print msg
    sys.exit(1)

def fixpath(path):
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))

def generate_csv_table(counts, csv_path):
    """
    Get the CSV file for the counts produced in get_all_counts()
    :param counts: list of count dicts
    :param csv_path: path to the csv file
    :return: generates file in csv_path
    """
    h1_ref = []
    h1_alt = []
    h2_ref = []
    h2_alt = []
    un_ref = []
    un_alt = []
    poss = []
    chroms = []
    ref = []
    alt = []
    for d in counts:
        h1_ref.append(d.get('h1', None).get('ref', None))
        h1_alt.append(d.get('h1', None).get('alt', None))
        h2_ref.append(d.get('h2', None).get('ref', None))
        h2_alt.append(d.get('h2', None).get('alt', None))
        un_ref.append(d.get('un', None).get('ref', None))
        un_alt.append(d.get('un', None).get('alt', None))
        poss.append(d.get('pos', None))
        chroms.append(d.get('chrom', None))
    df = pandas.DataFrame({
        'pos': poss,
        'chrom': chroms,
        'h1_ref': h1_ref,
        'h1_alt': h1_alt,
        'h2_ref': h2_ref,
        'h2_alt': h2_alt,
        'un_ref': un_ref,
        'un_alt': un_alt,
        'ref': ref,
        'alt': alt
    })
    df.to_csv(csv_path, index=False)

def get_all_counts(vcf_path, bam_path, ref):
    """
    get the count dict for all
    :param px_vcf:
    :param bam:
    :param fa:
    :return: a dict with all counts
    """

    fa = pyfasta.Fasta(ref)
    px_vcf = vcf.Reader(open(vcf_path))
    px_bam = pysam.Samfile(bam_path)

    counts = []
    for index, v in enumerate(px_vcf):
        record = get_counts_for_record(v, px_bam, fa)
        counts.append(record)
        if index % 20 == 0:
            print("Running record ", index + 1, " Record: ", record)
    return counts

def get_counts_for_record(vcf_rec, bam, fa):
    """
    Get the table h1/h2 ref/alt counts for a variant in vcf_record in the given bam file.
    :param vcf_rec:
    :param bam:
    :param fa:
    :return:
    """
    # 'ref': ref
    alleles = tk_io.get_record_alt_alleles(vcf_rec)
    ref = tk_io.get_record_ref(vcf_rec)
    r = get_allele_read_info(vcf_rec.CHROM, vcf_rec.POS, ref, alleles, 30, bam, fa)
    r['ref'] = ref
    r['alt'] = alleles[0]
    return r


def get_allele_read_info(chrom, pos, ref, alt_alleles, min_mapq, bam, 
                         reference_pyfasta, max_reads=2000, match = 1, 
                         mismatch = -4, gap_open = -6, gap_extend = -1):
    """
    The counts of each allele on each haplotype
    :param chrom:
    :param pos:
    :param ref:
    :param alt_alleles:
    :param min_mapq:
    :param bam:
    :param reference_pyfasta:
    :param max_reads:
    :param match:
    :param mismatch:
    :param gap_open:
    :param gap_extend:
    :return:
    """

    if not chrom.startswith('chr'):
        chrom = 'chr' + chrom

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
                                                              gap_open = gap_open,
                                                              gap_extend = gap_extend)
        for (allele_index, allele) in enumerate(all_alleles):
            if allele_index == allele_index_in_read:
                if read.mapq >= min_mapq:
                    counts[hap][allele_index] += 1

    counts_reformat = {k: {'ref': v[0], 'alt': v[1]} for (k,v) in counts.items()}
    counts_reformat['chrom'] = chrom
    counts_reformat['pos'] = pos
    return counts_reformat


__doc__ = '''
Get the counts of all - alt/ref, chrom, pos

Usage:
    count.py <ref_path> <vcf_path> <bam_path> <output_csv_path>

Arguments:
    ref_path            Path to the reference genome fasta
    vcf_path            Path to the VCF file for position reads
    bam_path            Path to the BAM file
    output_csv_path     Path of the CSV file you want for output

Options:
    -h --help   Show this message.
'''


def run():
    args = docopt.docopt(__doc__, version=" ")  # do not remove space
    ref_path = fixpath(args["<ref_path>"])
    vcf_path = fixpath(args["<vcf_path>"])
    bam_path = fixpath(args["<bam_path>"])
    output_csv_path = fixpath(args["<output_csv_path>"])

    counts = get_all_counts(vcf_path, bam_path, ref_path)
    generate_csv_table(counts, output_csv_path)


if __name__ == '__main__':
    run()
