echo "test start. "
if [ ! -d test_out ]; then mkdir test_out; fi

if [ ! -d test_out/marker_gene ]; then mkdir test_out/marker_gene; fi

echo "extrating marker genes..."
for i in EP1 FJAT15249.F50 GMI1000 SL2330; do 
    python bacphy/extract_ribosomal_protein_sequence.py -g tests/data/${i}.fna -a tests/data/${i}.gff -o test_out/marker_gene/${i}.fa; 
done

echo "organize sequence by genes..."
python bacphy/organize_by_gene.py -i test_out/marker_gene -o test_out/afa

echo "align genes by mafft..."
for i in  `ls test_out/afa | sed 's/.fa//'`; do echo "mafft --quiet test_out/afa/${i}.fa >test_out/afa/${i}.afa"; done |parallel

if [ ! -d test_out/concat ]; then mkdir test_out/concat; fi

python bacphy/concatenate.py -i test_out/afa -o test_out/concat/concated.afa -p test_out/concat/partition.txt

echo "building tree by concatenated sequence..."
if [ ! -d test_out/concat/iqtree ]; then mkdir test_out/concat/iqtree; fi
iqtree -quiet -s test_out/concat/concated.afa -spp test_out/concat/partition.txt -pre test_out/concat/iqtree/concat -m MFP
