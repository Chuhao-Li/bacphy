#!/usr/bin/env python3
__version__ = 'v0.1'
usage = '''organize sequence by gene
'''

import os
import argparse
from Bio import SeqIO

def a_gene_a_fasta(indir, outdir, seqtype="nt"):
    '''from a genome a fasta to a gene a fasta

    Parameters
    ------
    indir: str
        input directory, with marker genes organize by genome
    outdir:str
        output directory, with marker genes organize by gene

    Returns
    ------
    None
    '''
    if seqtype == "nt": 
        suffix = "fa" 
    elif seqtype == "aa":
        suffix = "faa"
    fastas = os.listdir(indir)
    fastas = [i for i in fastas if i.endswith(suffix)]
    collector = {}
    for fa in fastas:
        g = '.'.join(fa.split('.')[:-1])
        for rec in SeqIO.parse(indir.rstrip('/') + '/' + fa, 'fasta'):
            if rec.id in collector:
                collector[rec.id].append(f'>{g}\n{rec.seq}\n')
            else:
                collector[rec.id] = [f'>{g}\n{rec.seq}\n']
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for k,v in collector.items():
        with open(outdir.rstrip('/') + '/' + k + f'.{suffix}', 'w') as f:
            f.write(''.join(v))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.description = usage
    parser.add_argument('-i', '--indir', help="input directory with multiple sequence files")
    parser.add_argument('-o', '--outdir', help="output directory")
    parser.add_argument('-v', '--version', action='store_true', help="print version")
    args = parser.parse_args()
    a_gene_a_fasta(args.indir, args.outdir)

