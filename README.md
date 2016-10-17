# HackSeq 2016 - Project Team 8

Project Team 8 for HackSeq 2016 in Vancouver BC, Canada. The repository contains two projects.


## Project A: Finding Somatic Variants With Phasing (`/somatic/`)

Calling somatic mutation from tumor tissues only is challenge not only because you do not have a control to facilitate filtering out germline variants but it is difficult to differentiate low frequency somatic mutation from sequence noise/errors. In this study, we investigate whether we can leverage phasing information from reads to help differentiate somatic variants from germline alterations and sequencing errors.

### Code

This repository codebase dependes on 10xGenomics' `longranger` toolset. [Download and install longranger](http://support.10xgenomics.com/genome-exome/software/pipelines/latest/installation). It also depends on the linked-reads data from 10xGenomics.


#### `count/count.py` 

Get the counts of alt/ref (hap1, hap2, unphased, chrom, pos,  from VCF file.

```bash
python count.py [--bed=<bed>] <ref_path> <vcf_path> <bam_path> <output_csv_path>
```

return value: Writes to disk a CSV file (given by `output_csv_path`) with columns : `alt,chrom,filter,h1_alt,h1_ref,h2_alt,h2_ref,in_bed,pos,ref,un_alt,un_ref`

#### `somatic_probability.py`

Test run the somatic test on phased allele count data.

```bash
python somatic_test <count_file> <result_file>
```


## Project B: Metagenome


### Code
