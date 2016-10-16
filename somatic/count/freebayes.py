from __future__ import print_function # Python 2.x
import subprocess


def execute(cmd):
    """
    Execute a shell command and print output
    :param cmd:
    :return:
    """
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