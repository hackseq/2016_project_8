#!/usr/bin/env python
#
#
import os
import sys
import subprocess
import docopt
import shutil
import numpy
import scipy.stats
import pandas

__doc__ = '''
Test run the somatic test on phased allele count data

Usage:
    somatic_test <count_file> <result_file>

Arguments:
    count_file  Path to phased count data
    result_file Path to output file

Options:
    -h --help   Show this message.
'''


def error(msg):
    print msg
    sys.exit(1)

def fixpath(path):
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))


def main():
    args = docopt.docopt(__doc__, version=" ") # do not remove space
    count_path = fixpath(args["<count_file>"])
    result_path = fixpath(args["<result_file>"])

    count_data = pandas.read_csv(count_path)

    results = []
    for (i, r) in count_data.iterrows():
        print("row {0}".format(i))
        table = { 'h1': { 'ref': r.h1_ref, 'alt': r.h1_alt }, 'h2': { 'ref': r.h2_ref, 'alt': r.h2_alt }}

        ml = max_likelihood(table)
        results.append(ml)

    df = pandas.DataFrame(results)

    final = pandas.concat([count_data, df], axis=1)
    final.to_csv(result_path)


test_table = { 'h1': { 'ref': 99, 'alt': 5}, 'h2': {'ref': 99, 'alt': 7}}


def g0_likelihood(table, error_rate=None):
    """
    g0 - scenario where neither of the haplotype has somatic or germline variants.

    :param table:
    :param error_rate: you can give the error_rate, if None it runs get_optimal_error_rate()
    :return:
    """

    if error_rate is None:
        error_rate = get_optimal_error_rate(table, (0,0))

    n_v1 = table['h1']['ref'] + table['h1']['alt']
    v1 = scipy.stats.binom.pmf(table['h1']['alt'], n_v1, error_rate)

    n_v2 = table['h2']['ref'] + table['h2']['alt']
    v2 = scipy.stats.binom.pmf(table['h2']['alt'], n_v2, error_rate)
    
    return v1 * v2


def get_optimal_error_rate(table, mode):
    ''' Calculates the optimal max-likelihood error rate.
        @param mode 2-element tuple representing the error base rates for haplotypes 1 and 2 respectively.
                    Haplotype order is as in table. Examples: (0,0) for g0, (1,0) or (0,1) for heterozygous,
                    (1,1) for homozygous. Use None for somatic variant, i.e. (0, None) for h2 somatic.
    '''

    keys = ['h1', 'h2']
    calculation_dict = {}

    for i in range(2):
        ref = table[keys[i]]['ref']
        alt = table[keys[i]]['alt']
        error_mode = mode[i]

        if error_mode is None:
            num = 0.0
            denom = 0.0
        elif error_mode == 0:
            # if error_mode is 0 (expecting pure wilde type), then considers alt as the error_rate
            num = float(alt)
            denom = float(ref + alt)
        elif error_mode == 1:
            # if error_mode is 1 (expecting pure mutant), then considers ref as the error_rate
            num = float(ref)
            denom = float(ref + alt)

        calculation_dict[keys[i]] = {'num': num, 'denom': denom}

    totalNum = sum([calculation_dict[haplo]['num'] for haplo in ('h1', 'h2')]) 
    totalDenom = sum([calculation_dict[haplo]['denom'] for haplo in ('h1', 'h2')]) 

    if totalDenom == 0:
        optimal_error = 0.1
    else:
        optimal_error = float(totalNum) / float(totalDenom)

    if optimal_error > 0.1:
        optimal_error = 0.1  # cap off the error_rate so we do not miss a somatic variant (say 0.3)

    return optimal_error
        

def get_higher_vaf_haplotype(table):
    if table['h1']['alt'] == 0 and table['h1']['ref'] == 0:
        return ('h2', 'h1')
    elif table['h2']['alt'] == 0 and table['h2']['ref'] == 0:
        return ('h1', 'h2')
    else:
        h1_vaf = float(table['h1']['alt']) / float(table['h1']['ref'] + table['h1']['alt'])
        h2_vaf = float(table['h2']['alt']) / float(table['h2']['ref'] + table['h2']['alt'])
    
    higher_vaf_haplotype = 'h1' if h1_vaf > h2_vaf else 'h2'
    lower_vaf_haplotype = 'h1' if higher_vaf_haplotype == 'h2' else 'h2'

    return (higher_vaf_haplotype, lower_vaf_haplotype)

    
def g1_likelihood(table, error_rate = None, homozygous=False):
    ''' Outputs the likelihood value for Case 1: Germline variant.
        @param homozygous Boolean value indicating whether the germline variant is homozygous
        @param error_rate If float, calculates based on specific error rate. If None, chooses max-likelihood error rate.
    '''

    (higher_vaf_haplotype, lower_vaf_haplotype) = get_higher_vaf_haplotype(table)

    if error_rate is None:
        if homozygous:  # both haplotypes are pure 'alt'
            error_mode = (1,1)
        else:
            error_mode = (1,0) if higher_vaf_haplotype == 'h1' else (0,1)

        error_rate = get_optimal_error_rate(table, error_mode)
    
    # Haplotype 1
    n_v1 = table[higher_vaf_haplotype]['ref'] + table[higher_vaf_haplotype]['alt']
    v1 = scipy.stats.binom.pmf(table[higher_vaf_haplotype]['alt'], n_v1, 1-error_rate)
    
    v2_prob = 1-error_rate if homozygous else error_rate
    
    # Haplotype 2
    n_v2 = table[lower_vaf_haplotype]['ref'] + table[lower_vaf_haplotype]['alt']
    v2 = scipy.stats.binom.pmf(table[lower_vaf_haplotype]['alt'], n_v2, v2_prob)
    
    return v1 * v2


def g2_likelihood(table, error_rate = None):
    ''' Outputs the likelihood value for Case 2: Somatic variant.
    '''

    likelihoods = [(g2_likelihood_helper(table, float(var_percent)/100, error_rate), var_percent) for var_percent in range(5, 95, 5)]
    
    for variant_percent in range(5, 95, 5):
        variant_frequency = float(variant_percent)/100
        likelihood_value = g2_likelihood_helper(table, variant_frequency, error_rate)
    
    return(max(likelihoods))


def g2_likelihood_helper(table, variant_frequency, error_rate = 0.001):
    ''' Helper function which calculates the likelihood values for Case 2: Somatic variant
    '''
    
    (higher_vaf_haplotype, lower_vaf_haplotype) = get_higher_vaf_haplotype(table)

    if error_rate is None:
        error_mode = (None, 0) if higher_vaf_haplotype == 'h1' else (0, None)
        error_rate = get_optimal_error_rate(table, error_mode)
    
    # Haplotype 1: higher VAF haplotype
    n_v1 = table[higher_vaf_haplotype]['ref'] + table[higher_vaf_haplotype]['alt']
    v1 = scipy.stats.binom.pmf(table[higher_vaf_haplotype]['alt'], n_v1, variant_frequency)
        
    # Haplotype 2: lower VAF haplotype
    n_v2 = table[lower_vaf_haplotype]['ref'] + table[lower_vaf_haplotype]['alt']
    v2 = scipy.stats.binom.pmf(table[lower_vaf_haplotype]['alt'], n_v2, error_rate)
    
    #print('error_rate: {0}, frequency: {1}'.format(error_rate, variant_frequency))
    #print('v1: {0}, v2: {1}, total: {2}'.format(v1, v2, v1*v2))
    
    return v1 * v2
    
    
def likelihood_per_error_rate(table, error_rate=None):
    prior_g0 = 1.0
    prior_g1 = 1.0/1000
    prior_g2 = 1.0/10000
    
    g0 = g0_likelihood(table, error_rate) * prior_g0
    g1_het = g1_likelihood(table, error_rate, homozygous=False) * prior_g1
    g1_hom = g1_likelihood(table, error_rate, homozygous=True) * prior_g1
    g2 = g2_likelihood(table, error_rate) * prior_g2
    
    return {"g0": (g0, 0.0), 
            "g1_het": (g1_het, 0.0),
            "g1_hom": (g1_hom, 0.0),
            "g2": g2 }
    

def max_likelihood(table):

    # for now, write NaN's in all cases with 0's
    if table['h1']['ref'] == 0 and table['h2']['ref'] == 0 and table['h1']['alt'] == 0 and table['h2']['alt'] == 0:
        cols = {'g0': 'NaN',
                'g1_het': 'NaN',
                'g1_hom': 'NaN',
                'g2': 'NaN',
                'best_model': 'NaN',
                'param': 'NaN'
                }
    else:
        likelihood_dict = likelihood_per_error_rate(table, error_rate=None)
        ((max_likelihood, param), model_name) = max((likelihood, name) for (name, likelihood) in likelihood_dict.items())
        
        z = sum(prob for (prob, param) in likelihood_dict.values())
        cols = {    'g0': likelihood_dict['g0'][0]/z, 
                    'g1_het': likelihood_dict['g1_het'][0]/z, 
                    'g1_hom': likelihood_dict['g1_hom'][0]/z, 
                    'g2': likelihood_dict['g2'][0]/z,
                    'best_model': model_name, 
                    'param': param }

    #print("Model: {0} had likelihood: {1}".format(model_name, max_likelihood))
    
    return cols


if __name__ == '__main__':
    main()

# scipy.stats.binom.pmf(3, 10, 0.1)
# print test_table
# max_likelihood(test_table)
