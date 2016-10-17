# Somatic Mutation from Separated Haplotypes (SMUSH)

## Team (in alphabetical order) 

Amanjeev Sethi, Eric Zhao, Hua Ling, Patrick Marks (Project lead), Peng Zhang, Samantha Kohli


## Scientific Abstract

Currently, accurate somatic mutation calling relies upon comparison of tumor to matched normal. The matched normal enables filtering of germline variants and differentiation of low frequency somatic mutations from sequence noise/errors. However, matched normal tissue is often unavailable. In this study, we investigate whether phasing information from linked reads can help differentiate somatic variants from germline alterations and sequencing errors. 

Whole exome sequencing of the HCC1954 BrCa cell line was performed at 0%, 25%, 50%, 75% and 100% tumor purity. Library construction was performed using 10x genomics linked reads. Samples were sequenced at mean coverage of 200x. Sequence alignment and variant calling were performed using longranger and FreeBayes. Read counts of reference and alternate alleles were obtained for each haplotype after excluding reads with mapping quality < 30 to calculate variant allele fractions (VAF). To differentiate between wild type, germline variants, and somatic mutations, model selection was performed using a maximum-likelihood estimate accounting for haplotype phasing, mutant read capture rates, and sequencing error. 

To benchmark the accuracy of model selection, 31 somatic variants validated by comparison of 0% and 100% tumor content samples were manually reviewed in IGV and selected as a testing subset. Employing the SMUSH model on this subset correctly identified all variants as somatic in the 25-100% samples and correctly identified 30/31 variants as wild-type in the 0% sample.


## Lay Abstract

Sequencing of tumor genomes has become a critical component of cancer care. In the clinical setting, the normal cell sample is usually absent due to logistics and/or monetary limitations, but it important for identifying mutations that are specific to the cancer sample. 10x Genomics has specialized genetic data which allows the assignment of mutations to a specific category (haploytype) within the genome by providing adding long-range information to short read sequencing. Cancer-causing variation is generally and partially associated with only one of these categories (haplotypes). The long-range information helps identify and separate out tumor-specific variation from germline and sequencing error, in the absence of normal cell data to compare it with. In this project, we attempt that analysis using an existing variant data for HCC1954 and given long-range data from 10x.
