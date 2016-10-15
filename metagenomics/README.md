# Metagenomics Project

We have done some basic prep of the metagenomic sequencing data to help get things rolling. Basic 10x barcode processing is performed to extract barcodes from reads & apply barcode error correction.The reads are then fed into a simple de Bruijn graph assembler that keeps track of the set of barcodes that contribute a kmer to each unbranched path in the graph.  

This process generated the following outputs:

* asm.gfa - a representation of the assembly graph in [GFA](https://github.com/GFA-spec/GFA-spec/blob/master/GFA-spec.md) format.  This file can be easily loaded into the [Bandage](https://rrwick.github.io/Bandage/) assembly viewer, to inspect the structure of the graph.

* stats.tsv - a simple TSV file with high-level statistics about each contig in the graph. Columns are:
  * contig_id
  * length of contig sequence
  * number of barcode on contig
  * number of outbound edges on to the left
  * number of outbound edges to the right
  * contig sequence

* node_bcs.mtx.gz - the set of barcodes contributing to each contig of the graph, encoded as a sparse adjacency matrix in [MatrixMarket](http://people.sc.fsu.edu/~jburkardt/data/mm/mm.html) format. This file can be easily read into a Scipy sparse ccoo matrix with the [mmread](http://docs.scipy.org/doc/scipy/reference/generated/scipy.io.mmread.html) command like so: 
```
In [11]: gr = scipy.io.mmread("node_bcs.mtx.gz")

In [12]: gr
Out[12]:
<107550x4792320 sparse matrix of type '<type 'numpy.float64'>'
        with 18229421 stored elements in COOrdinate format>
```

Bandage screenshot:
![alt text](https://github.com/10XDev/hackseq_dna/raw/master/metagenomics/banadge-example.png "Logo Title Text 1")
