# bacphy v0.2

### Introduction
Bacphy(BACterial PHYlogenetic analysis tools) aims to build a robust phylogenetic tree by bacterial genomes on different taxonomic level. 
According to the taxonomic level of our input bacterial genomes, we can use different sequence to build a tree. For bacterial from different families, 16S rRNA may be sufficient. While in strain level, we may need more sequence to ensure enough phylogenetic information. Here are some options:

1. small sub-unit rRNA(16S rRNA) sequence
2. animo acid sequence of marker genes
3. nucleotide sequence of marker genes
4. all protein conding sequences of core genome
5. whole genome(including intergenic region)

Here, 120 single-copy marker genes are used . This pipeline includes the following steps: 

1. extract sequence from genome sequence. 
2. align sequence. 
3. build tree. 

### Installation: 
This program is developed in Ubuntu 20.04 and has not been tested in other environment. 

pre-requests: 
- python v3.7.3
- python package pyfaidx, biopython
- Prodigal V2.6.3
- HMMER 3.3
- mafft v7.455
- iqtree v1.6.12

to install the main program: 
``` bash
git clone https://github.com/Chuhao-Li/bacphy.git

cd bacphy

python3 setup.py install
```

### Test the program

``` bash
bacphy -g tests/data/genome/ -o test_out 
```

### more infomation:

New feature of v0.2 
- Easier installing and running
- More marker genes
- No need of annotation file 

To do list: 
- Raise warning or error if the sequence quality is low. 
- Raise warning if multi-copy of marker genes exists

