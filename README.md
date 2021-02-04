# bacphy v0.1

### Introduction
Bacphy(BACterial PHYlogenetic analyais tools) aims to build a robust phylogenetic tree by bacterial genomes on different taxonomic level. 
According to the taxonomic level of our input bacterial genomes, we can use different sequence to build a tree. For bacterial from different families, 16S rRNA may be sufficient. While in strain level, we may need more sequence to ensure enough phylogenetic information. Here are some choices:

1. small sub-unit rRNA(16S rRNA) sequence
2. animo acid sequence of marker genes
3. nucleotide sequence of marker genes
4. all protein conding sequences of core genome
5. whole genome(including intergenic region)

In many cases, nucleotide sequence of ribosomal protein genes is enough to differentiate different species. So, in the v0.1 of bacphy, I will implement a pipeline to build trees by nucleotide sequence of marker genes. This pipeline includes the following steps: 

1. extract sequence from genome sequence. 
2. align sequence. 
3. build tree. 

### Installation: 
This program is developed in Ubuntu 20.04 and has not been tested in other environment. 

pre-requests: 
- python v3.7.3
- python package pyfaidx, biopython
- parallel 20160622
- mafft v7.455
- iqtree v1.6.12

to install the main program: 
``` bash
git clone https://github.com/Chuhao-Li/bacphy.git
```

### Test the program
If the program is run sucessfully, a directory named test_out will be generated and its content should be the same as tests/test1_out

``` bash
cd bacphy

bash tests/script/test1.sh
```

### Cookbook
Here, we descript what happened in the test procedure. 

Input files are genome sequence file and annotation file of 4 Ralstonia solanacearum strains: 
```
tests/data/
├── EP1.fna
├── EP1.gff
├── FJAT15249.F50.fna
├── FJAT15249.F50.gff
├── GMI1000.fna
├── GMI1000.gff
├── SL2330.fna
└── SL2330.gff
```

Create output directory. 
``` bash
echo "test start. "
if [ ! -d test_out ]; then mkdir test_out; fi
if [ ! -d test_out/marker_gene ]; then mkdir test_out/marker_gene; fi
```

Extract ribosomal protein genes. 
``` bash
echo "extrating marker genes..."
for i in EP1 FJAT15249.F50 GMI1000 SL2330; do
    python bacphy/extract_ribosomal_protein_sequence.py -g tests/data/${i}.fna -a tests/data/${i}.gff -o test_out/marker_gene/${i}.fa;
done

echo "organize sequence by genes..."
python bacphy/organize_by_gene.py -i test_out/marker_gene -o test_out/afa
```

Align sequences. 
``` bash
echo "align genes by mafft..."
for i in  `ls test_out/afa | sed 's/.fa//'`; do echo "mafft --quiet test_out/afa/${i}.fa >test_out/afa/${i}.afa"; done |parallel
```

Concatenate all aligned sequence of the same genome, and build tree. 
``` bash
if [ ! -d test_out/concat ]; then mkdir test_out/concat; fi

python bacphy/concatenate.py -i test_out/afa -o test_out/concat/concated.afa -p test_out/concat/partition.txt

echo "building tree by concatenated sequence..."
if [ ! -d test_out/concat/iqtree ]; then mkdir test_out/concat/iqtree; fi
iqtree -quiet -s test_out/concat/concated.afa -spp test_out/concat/partition.txt -pre test_out/concat/iqtree/concat -m MFP
```
