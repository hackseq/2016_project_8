#!/usr/bin/env python
#
#
"""
Get the counts of alt/ref (hap1, hap2, unphased, chrom, pos,  from VCF file

Usage:
    count.py [--bed=<bed>] <ref_path> <vcf_path> <bam_path> <output_csv_path>

Arguments:
    bed                 Optional bed file for filtering
    ref_path            Path to the reference genome fasta
    vcf_path            Path to the VCF file for position reads
    bam_path            Path to the BAM file
    output_csv_path     Path of the CSV file you want for output

Options:
    -h --help   Show this message.

"""
import os
import sys

import docopt
import pandas

import pyfasta
import pysam
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import vcf


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
    data_frame = pandas.DataFrame(counts)
    data_frame.to_csv(csv_path, index=False)

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
    alleles = tk_io.get_record_alt_alleles(vcf_rec)
    ref = tk_io.get_record_ref(vcf_rec)
    r = get_allele_read_info(vcf_rec.CHROM, vcf_rec.POS, ref, alleles, 30, bam, fa)
    r['ref'] = ref
    r['alt'] = alleles[0]
    filter_col_vcf = ':'.join(vcf_rec.FILTER)
    if filter_col_vcf is None or filter_col_vcf.strip() == "":
        filter_col_vcf = "PASS"  # PASS does not return a list and gets dropped
    r['filter'] = filter_col_vcf
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

    counts_reformat = {}
    for k, v in counts.items():
        counts_reformat[k + '_' + 'ref']  = v[0]
        counts_reformat[k + '_' + 'alt'] = v[1]
    counts_reformat['chrom'] = chrom
    counts_reformat['pos'] = pos
    return counts_reformat

def annotate_bed_info(counts, bed_file):
    regs = tk_io.get_target_regions(open(bed_file))

    for c in counts:
        in_bed = False
        if regs.has_key(c['chrom']):
            in_bed = regs[c['chrom']].contains_point(c['pos'])

        c['in_bed'] = in_bed


__doc__ = '''
Get the counts of all - alt/ref, chrom, pos

Usage:
    count.py [--bed=<bed>] <ref_path> <vcf_path> <bam_path> <output_csv_path>

Arguments:
    bed                 Optional bed file for filtering
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

    bed_path = None
    if args["--bed"] is not None:
        bed_path = fixpath(args["--bed"])

    counts = get_all_counts(vcf_path, bam_path, ref_path)

    if bed_path is not None:
        annotate_bed_info(counts, bed_path)

    generate_csv_table(counts, output_csv_path)


if __name__ == '__main__':
    run()
