#!/usr/bin/env python3

usage = '''*.afa -> concated.afa
1. {genome1: [gene1], genome2: [gene1], ...}
2. {genome1: [gene1, gene2], genome2: [gene1, gene2], ...}
...

Args: 
    *.afa

Returns:
    concated.afa
'''

import os
import argparse
from Bio import SeqIO

def concatenate(indir, outdir, suffix="afa"): # suffix: output suffix
    outdir = outdir.rstrip('/')
    indir = indir.rstrip('/')
    outfa = f'{outdir}/concatenated.{suffix}'
    out_par = f'{outdir}/partition.txt'
    fastas = [i for i in os.listdir(indir) if i.endswith(f'.{suffix}')]
    
    out = {}
    
    fa = fastas[0]
    gene = fa.replace(f'.{suffix}', '')
    partition = []
    for rec in SeqIO.parse(f'{indir}/{fa}', 'fasta'):
        out[rec.id] = [str(rec.seq)]
    else:
        partition.append((gene, len(rec.seq)))
    
    for fa in fastas[1:]:
        gene = fa.replace(f'.{suffix}', '')
        for rec in SeqIO.parse(f'{indir}/{fa}', 'fasta'):
            out[rec.id].append(str(rec.seq))
        else:
            partition.append((gene, len(rec.seq)))
    
    # write concated sequence
    with open(outfa, 'w') as f:
        for k,v in out.items():
            concated = ''.join(v)
            f.write(f'>{k}\n{concated}\n')
    
    # write partition
    model = "protein" if suffix == "afaa" else "DNA"
    # if use -TEST, iqtree will ignore the "protein" model. see https://github.com/Cibiv/IQ-TREE/issues/153
    with open(out_par, 'w') as f:
        start = 1
        for gene, l in partition:
            end = str(start+l-1)
            f.write(f'{model}, {gene} = {start}-{end}\n')
            start += l

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.description = usage
    parser.add_argument('-i', '--indir', help="input directory with multiple sequence files")
    parser.add_argument('-o', '--outdir', help="output aligned fasta file")
    parser.add_argument('-v', '--version', action='store_true', help="print version")
    args = parser.parse_args()
    concatenate(args.indir, args.outdir)
