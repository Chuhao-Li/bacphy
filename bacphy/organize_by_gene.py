#!/usr/bin/env python3
__version__ = 'v0.1'
usage = '''organize sequence by gene
'''

import os
import argparse
from Bio import SeqIO

def main(args):
    fastas = os.listdir(args.indir)
    fastas = [i for i in fastas if i.split('.')[-1] in {'fa', 'fna', 'fasta'}]
    collector = {}
    for fa in fastas:
        g = '.'.join(fa.split('.')[:-1])
        for rec in SeqIO.parse(args.indir.rstrip('/') + '/' + fa, 'fasta'):
            if rec.id in collector:
                collector[rec.id].append(f'>{g}\n{rec.seq}\n')
            else:
                collector[rec.id] = [f'>{g}\n{rec.seq}\n']
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    for k,v in collector.items():
        with open(args.outdir.rstrip('/') + '/' + k + '.fa', 'w') as f:
            f.write(''.join(v))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.description = usage
    parser.add_argument('-i', '--indir', help="input directory with multiple sequence files")
    parser.add_argument('-o', '--outdir', help="output directory")
    parser.add_argument('-v', '--version', action='store_true', help="print version")
    args = parser.parse_args()
    main(args)

