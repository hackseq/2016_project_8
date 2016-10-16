#!/usr/bin/env python
#
#
"""
Main entry point
"""
import docopt
from freebayes import run_freebayes
from count import get_all_counts, generate_csv_table

__doc__ = '''
Get the counts of all - alt/ref, chrom, pos

Usage:
    main.py <ref_path> <vcf_path> <bam_path> <output_csv_path>

Arguments:
    ref_path            Path to the reference genome fasta
    vcf_path            Path to the VCF file for position reads
    bam_path            Path to the BAM file
    output_csv_path     Path of the CSV file you want for output

Options:
    -h --help   Show this message.
'''

def error(msg):
    print msg
    sys.exit(1)

def fixpath(path):
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))

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