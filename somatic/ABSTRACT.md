# Somatic Mutation from Separated Haplotypes (SMUSH)

## Team (in alphabetical order) 

Amanjeev Sethi, Eric Zhao, Hua Ling, Patrick Marks (Project lead), Peng Zhang, Samantha Kohli


## Scientific Abstract

Currently, accurate somatic mutation calling relies upon comparison of tumor to matched normal. The matched normal enables filtering of germline variants and differentiation of low frequency somatic mutations from sequence noise/errors. However, matched normal tissue is often unavailable. In this study, we investigate whether phasing information from linked reads can help differentiate somatic variants from germline alterations and sequencing errors. 

Whole exome sequencing of the HCC1954 BrCa cell line was performed at 0%, 25%, 50%, 75% and 100% tumor purity. Library construction was performed using 10x genomics linked reads. Samples were sequenced at mean coverage of 200x. Sequence alignment and variant calling were performed using longranger and FreeBayes. Read counts of reference and alternate alleles were obtained for each haplotype after excluding reads with mapping quality < 30 to calculate variant allele fractions (VAF).

To differentiate between wild type, germline variants, and somatic mutations, we devised the SMUSH algorithm. Haplotype-specific VAFs were modeled as binomial random processes. model selection was performed using a maximum-likelihood estimate accounting for haplotype phasing, mutant read capture rates, and sequencing error. 

To benchmark the accuracy of model selection, 31 somatic variants validated by comparison of 0% and 100% tumor content samples were manually reviewed in IGV and selected as a testing subset. Employing the SMUSH model on this subset correctly identified all variants as somatic in the 25-100% samples and correctly identified 30/31 variants as wild-type in the 0% sample.



## Lay Abstract

Sequencing of tumor genomes has become a critical component of cancer care. Sequencing a non-tumor sample from the patient is  important for identifying mutations that are specific to the cancer sample, but in the clinical setting, the normal sample is often not available due to logistics and/or monetary limitations. 10x Genomics technology adds long-range information to the short read sequencing used to interrogate cancer genomes. Long-range information enables the assignment of mutations to a specific category (haploytype) within the genome. Mutations specific to the cancer will show a distintive pattern across haplotypes, which can be used to separate them from non-cancer mutations, without separately sequencing the normal sample. In this project, we build a prototype method and apply it to a standard cancer test sample HCC1954 sequenced with 10x Genomics Chromium technology.
