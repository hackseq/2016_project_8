# Somatic Mutation from Separated Haplotypes (SMUSH)

## Team (in alphabetical order) 

Amanjeev Sethi, Eric Zhao, Hua Ling, Patrick Marks (Project lead), Peng Zhang, Samantha Kohli


## Scientific Abstract

Calling somatic mutation from tumor tissues only is challenge not only because you do not have a control to facilitate filtering out germline variants but it is difficult to differentiate low frequency somatic mutation from sequence noise/errors. In this study, we investigate whether we can leverage phasing information from reads to help differentiate somatic variants from germline alterations and sequencing errors. 

HCC1954 BrCa cell lines were WES/WGS sequenced using ilmn *** at 0, 25, 50, 75 and 100% tumor purity. 10x genomics linked reads was used for library preparation (Pat). Samples were sequenced at mean coverage of **x. Freebayes and longranger ** were used for alignment and **variant calling? Read counts of reference and alternate alleles were obtained for each haplotype after excluding reads with MAQP<30 to calculate variant allele fractions (VAF). Likelihood of being noise, germline and somatic variants were calculated based on binomial distribution. …. The state of variants are determined using maximum likelihood model (Eric). (shall we say something about looking at VAF in different tiers of the cellline?)Variants with high VAF on one haplotype and low VAF on the other are manually reviewed in IGV to optimize our somatic variant calling algorithm. A quick sanity check was applied (compare the proportion of VAF from two haplotypes).


## Lay Abstract

In most studies, the normal cell sample is usually absent. This could be due to logistics or monetary limitations. 10x Genomics has a specialized genetic data which allows the categorization within the genome by providing long-range information. This long-range information helps identify and separate out the clusters (which are categorized) and helps with the analysis of tumor specific genetic variation. Ultimately, in the absence of normal cell data to compare with, we can potentially make use of the 10x's long-range data. In this project, we attempt that using an existing variant data for HCC1954. 