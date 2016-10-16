
# coding: utf-8

# In[52]:

import numpy
import scipy.stats


# In[266]:


test_table = { 'h1': { 'ref': 99, 'alt': 5}, 'h2': {'ref': 99, 'alt': 7}}


# In[267]:

def g0_likelihood(table, error_rate=0.001):
    n_v1 = table['h1']['ref'] + table['h1']['alt']
    v1 = scipy.stats.binom.pmf(table['h1']['alt'], n_v1, error_rate)

    n_v2 = table['h2']['ref'] + table['h2']['alt']
    v2 = scipy.stats.binom.pmf(table['h2']['alt'], n_v2, error_rate)
    
    return v1 * v2
    
    
def g1_likelihood(table, error_rate = 0.001, homozygous=False):
    ''' Outputs the likelihood value for Case 1: Germline variant.
        @param homozygous Boolean value indicating whether the germline variant is homozygous
    '''
    
    h1_vaf = float(table['h1']['alt']) / float(table['h1']['ref'] + table['h1']['alt'])
    h2_vaf = float(table['h2']['alt']) / float(table['h2']['ref'] + table['h2']['alt'])
    
    higher_vaf_haplotype = 'h1' if h1_vaf > h2_vaf else 'h2'
    lower_vaf_haplotype = 'h1' if higher_vaf_haplotype == 'h2' else 'h2'
    
    # Haplotype 1
    n_v1 = table[higher_vaf_haplotype]['ref'] + table[higher_vaf_haplotype]['alt']
    v1 = scipy.stats.binom.pmf(table[higher_vaf_haplotype]['alt'], n_v1, 1-error_rate)
    
    v2_prob = 1-error_rate if homozygous else error_rate
    
    # Haplotype 2
    n_v2 = table[lower_vaf_haplotype]['ref'] + table[lower_vaf_haplotype]['alt']
    v2 = scipy.stats.binom.pmf(table[lower_vaf_haplotype]['alt'], n_v2, v2_prob)
    
    return v1 * v2


def g2_likelihood(table, error_rate = 0.001):
    ''' Outputs the likelihood value for Case 2: Somatic variant.
    '''
    
    likelihoods = [g2_likelihood_helper(table, float(var_percent)/100, error_rate) for var_percent in range(5, 95, 5)]
    
    for variant_percent in range(5, 95, 5):
        variant_frequency = float(variant_percent)/100
        likelihood_value = g2_likelihood_helper(table, variant_frequency, error_rate)
    
    return(max(likelihoods))


def g2_likelihood_helper(table, variant_frequency, error_rate = 0.001):
    ''' Helper function which calculates the likelihood values for Case 2: Somatic variant
    '''
    
    h1_vaf = float(table['h1']['alt']) / float(table['h1']['ref'] + table['h1']['alt'])
    h2_vaf = float(table['h2']['alt']) / float(table['h2']['ref'] + table['h2']['alt'])
    
    higher_vaf_haplotype = 'h1' if h1_vaf > h2_vaf else 'h2'
    lower_vaf_haplotype = 'h1' if higher_vaf_haplotype == 'h2' else 'h2'
    
    # Haplotype 1: higher VAF haplotype
    n_v1 = table[higher_vaf_haplotype]['ref'] + table[higher_vaf_haplotype]['alt']
    v1 = scipy.stats.binom.pmf(table[higher_vaf_haplotype]['alt'], n_v1, variant_frequency)
        
    # Haplotype 2: lower VAF haplotype
    n_v2 = table[lower_vaf_haplotype]['ref'] + table[lower_vaf_haplotype]['alt']
    v2 = scipy.stats.binom.pmf(table[lower_vaf_haplotype]['alt'], n_v2, error_rate)
    
    #print('error_rate: {0}, frequency: {1}'.format(error_rate, variant_frequency))
    #print('v1: {0}, v2: {1}, total: {2}'.format(v1, v2, v1*v2))
    
    return v1 * v2
    
    
def likelihood_per_error_rate(table, error_rate=0.001):
    prior_g0 = 1
    prior_g1 = 1
    prior_g2 = 1
    
    g0 = g0_likelihood(table, error_rate) * prior_g0
    g1_het = g1_likelihood(table, error_rate, homozygous=False) * prior_g1
    g1_hom = g1_likelihood(table, error_rate, homozygous=True) * prior_g1
    g2 = g2_likelihood(table, error_rate) * prior_g2
    
    return {"g0": g0, 
            "g1_het": g1_het,
            "g1_hom": g1_hom,
            "g2": g2}
    

def max_likelihood(table):
    maximum = -1
    maximum_model = None
    maximum_likelihood_error = -1
    
    for i in range(1, 100):
        error_rate = float(i)/1000
        likelihood_dict = likelihood_per_error_rate(table, error_rate)
        
        max_likelihood = max(likelihood_dict.values())
        max_model = [model_name for model_name in likelihood_dict if likelihood_dict[model_name] == max_likelihood][0]
        if max_likelihood > maximum:
            maximum = max_likelihood
            maximum_model = max_model
            maximum_likelihood_error = error_rate
            
    print("Model: {0}, with error rate: {1}, had likelihood: {2}".format(maximum_model, maximum_likelihood_error, maximum))


# In[272]:

scipy.stats.binom.pmf(3, 10, 0.1)


# In[269]:

print test_table
max_likelihood(test_table)


# In[ ]:




# In[ ]:



