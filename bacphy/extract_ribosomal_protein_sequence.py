#!/usr/bin/env python3

version = '0.1'
usage = f'''
extract ribosomal protein sequence from fasta
version: {version}

usage: 
python bacphy/extract_ribosomal_protein_sequence.py \\
        -g test/data/Ralstonia_solanacearum_EP1_genomic.fna \\
        -a test/data/Ralstonia_solanacearum_EP1_genomic.gff \\
        -o test_out

input: a genome sequence file(*.fasta) and corresponding annotation file(*.gff)
output: extracted ribosomal protein sequence(*.fasta)
'''

import re
import argparse
from pyfaidx import Fasta
from Bio import SeqIO

def check_single_copy(rp_dict_count):
    '''ribosomal protein copy number statistics. 

    Args: 
        rp_dict_count: dict, {'rp1': [], 'rp2': [], ...}
    Returns: 
        dict, {'loss': ['rp1', ...], 
               'single': ['rp2', ...], 
               'duplicated': ['rp10', ...], 
               'present':[...]
              }
    '''
    rp_set = set(rp_dict_count.values())
    a = {'loss': [], 'single': [], 'duplicated': [], 'present':[]}
    for k,v in rp_dict_count.items():
        if v == 0:
            a['loss'].append(k)
        elif v== 1:
            a['single'].append(k)
            a['present'].append(k)
        else:
            a['duplicated'].append(k)
            a['present'].append(k)
    return a

def get_ribosomal_coord(gff_file, rp_list):
    '''get coordinate of ribosomal protein gene.

    Args: 
        gff_file: string, path of gff file
        rp_list: list, a list of ribosomal protein keywords
    Returns: 
        dict, {'rp1':gffline1, 'rp2':gffline2, ...}
    '''
    rp_dict = dict.fromkeys(rp_list, '')
    rp_dict_count = dict.fromkeys(rp_list, 0)
    for line in open(gff_file):
        # product = re.search(r'(?<=product=)[^;\n]*', line)
        product = re.search(r'(?<=gene=)[^;\n]*', line)
        if product:
            product = product.group()
            if product in rp_dict:
                rp_dict[product] = line
                # print(line)
                  # in this way, if a rp is not single copy, 
                  # only the last one will be extracted. 
                rp_dict_count[product] += 1
        else:
            continue
    # a = check_single_copy(rp_dict_count)
    return rp_dict

def get_fasta_by_coord(line, fasta):
    '''subset sequence from genome sequence file. 

    Args: 
        line: string, a gff line
        fasta: string, path of genome sequence file
    Returns: 
        string, the ribosomal protein gene sequence
    '''
    if not line:
        return ''
    genes = Fasta(fasta)
    line = line.split('\t')
    chro = line[0]
    startPos = int(line[3])
    endPos = int(line[4])
    strand = line[6]
    sequence = genes[chro][startPos-1:endPos]
    if strand == '-':
        sequence = sequence.reverse.complement
    return sequence.seq

def main(args):
    if args.version:
        print(version)
        exit()
    #rp_list = ["30S ribosomal protein S" + str(i) for i in range(1, 22)] + \
    #          ["50S ribosomal protein L" + str(i) for i in range(1, 37) if i not in [7, 8, 12, 26]] + \
    #          ['50S ribosomal protein L7/L12'] 
    rp_list = open('data/single_copy_genes.list').read().strip().splitlines()
    rp_list = [i.split('\t')[1] for i in rp_list]
    coords = get_ribosomal_coord(args.anno, rp_list)
    with open(args.outfile, 'w') as f:
        for rp, coord in coords.items():
            seq = get_fasta_by_coord(coord, args.genome)
            f.write(f'>{rp}\n{seq}\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.description = usage
    parser.add_argument('-g', '--genome', help="genome sequence file(*.fasta or *.fna or *.fa)")
    parser.add_argument('-a', '--anno', help="genome annotation file(*.gff)")
    parser.add_argument('-o', '--outfile', help="output file")
    parser.add_argument('-v', '--version', action='store_true', help="print version")
    args = parser.parse_args()
    main(args)
