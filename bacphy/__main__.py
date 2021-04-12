#!/usr/bin/env python3
import os
import sys
import argparse
from multiprocess import Pool
from tqdm import tqdm
from Bio import SeqIO
from bacphy.organize_by_gene import a_gene_a_fasta
from bacphy.concatenate import concatenate

def parse_args():
    usage = 'bacterial phylogenetic analysis tools'
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.description = usage
    parser.add_argument('-g', '--genome_dir', help="directory including genome sequence file(*.fna)")
    parser.add_argument('-t', '--threads', default=8, type=int, help="num of processes to use")
    parser.add_argument('-f', '--force', action="store_true", help="force to over-write existed prodigal and hmmsearch result")
    parser.add_argument('-o', '--outdir', default='bacphy_out', help="output directory")
    args = parser.parse_args()
    return args

def extract_marker_gene(indir, outdir, tmp, data_dir, threads, force):
    '''extract 120 marker genes from genome sequence

    Parameters
    ------
    seqs: str
        directory with genome sequence

    Returns
    ------
    dict
        {
        genome: {
            marker gene name: marker gene sequence, ...}
            , ...}
    '''
    genomes = os.listdir(indir)
    prodigal_tmp = f'{tmp}/prodigal'
    if not os.path.exists(prodigal_tmp):
        os.mkdir(prodigal_tmp)
    commands = []
    for i in genomes:
        prefix = i.replace('.fna', '')
        command = ['prodigal']
        command.extend(['-a', f'{prodigal_tmp}/{prefix}.faa'])
        command.extend(['-d', f'{prodigal_tmp}/{prefix}.cds'])
        command.extend(['-o', f'{prodigal_tmp}/{prefix}.gff'])
        command.extend(['-c', '-q', '-f', 'gff', '-i', f'{indir}/{i}'])
        if force or not os.path.exists(f'{prodigal_tmp}/{prefix}.faa'):
            commands.append(' '.join(command))
    print('running prodigal...')
    with Pool(processes=threads) as p:
        with tqdm(total=len(commands)) as pbar:
            for i in p.imap(os.system, commands):
                pbar.update()
    # hmmsearch
    hmmsearch_tmp = f'{tmp}/hmmsearch'
    if not os.path.exists(hmmsearch_tmp):
        os.mkdir(hmmsearch_tmp)
    commands = []
    for g in genomes:
        prefix = g.replace('.fna', '')
        for i in ['pfam', 'tigr']:
            command = ['hmmsearch', '-o', f'{hmmsearch_tmp}/{prefix}.txt']
            command.extend(['--tblout', f'{hmmsearch_tmp}/{prefix}.{i}.tsv'])
            command.extend(['--cut_nc'])
            command.extend([f'{data_dir}/hmm_db/bac120.{i}.hmm'])
            command.extend([f'{prodigal_tmp}/{prefix}.faa'])
            if force or not os.path.exists(f'{hmmsearch_tmp}/{prefix}.{i}.tsv'):
                commands.append(' '.join(command))
    print('running hmmsearch...')
    with Pool(processes=threads) as p:
        with tqdm(total=len(commands)) as pbar:
            for i in p.imap(os.system, commands):
                pbar.update()

    # extract single copy gene
    marker_genes = {}
    for g in genomes:
        prefix = g.replace('.fna', '')
        prot_to_hmm = {}
        for i in ['pfam', 'tigr']:
            db = f'{hmmsearch_tmp}/{prefix}.{i}.tsv'
            for line in open(db):
                if line.startswith('#'):
                    continue
                sline = line.split() # 这个函数相当于re.split(r'\s+', line)
                prot, hmm, score = sline[0], sline[3], float(sline[5])
                if prot not in prot_to_hmm:
                    prot_to_hmm[prot] = (hmm, score)
                elif score >prot_to_hmm[prot][1]:
                    prot_to_hmm[prot] = (hmm, score)
        hmm_to_prot = {} # {hmm: [prot1, prot2...]}
        for prot,hmm in prot_to_hmm.items():
            hmm = hmm[0]
            if hmm not in hmm_to_prot:
                hmm_to_prot[hmm] = [prot]
            else:
                hmm_to_prot[hmm].append(prot)
        marker_genes[prefix] = hmm_to_prot
    
    hmm_accs = [i.split('\t')[0] for i in open(f'{data_dir}/gtdbtk.bac120.marker_info.tsv').read().splitlines()[1:]]
    # 输出marker genes统计表
    with open(f'{outdir}/stat.tsv', 'w') as f:
        print('genome\tsingle\tmulti\tabsent', file=f)
        for g, v in marker_genes.items():
            single = 0
            multi = 0
            absent = 0
            for hmm in hmm_accs:
                if hmm in v:
                    if len(v[hmm]) > 1:
                        multi += 1
                    elif len(v[hmm]) == 1:
                        single += 1
                else:
                    absent += 1
            print(f'{g}\t{single}\t{multi}\t{absent}', file=f)

    # 输出单拷贝基因序列
    out_dict = {}
    for g, v in marker_genes.items():
        out_dict[g] = {}
        out = open(f'{outdir}/{g}.fa', 'w')
        cds = {rec.id:rec for rec in SeqIO.parse(f'{prodigal_tmp}/{g}.cds', 'fasta')}
        for acc in hmm_accs:
            seq = ''
            if acc in v:
                v1 = v[acc]
                seq = str(cds[v1[0]].seq)
            outfa = f'>{acc}\n{seq}\n'
            out.write(outfa)
            out_dict[g][acc] = outfa
    return out_dict

def align(dir2, dir3, threads):
    '''align all fa files in dir2'''

    commands = []
    for i in os.listdir(dir2):
        if not i.endswith('.fa'):
            continue
        prefix = i.replace('.fa', '')
        _in = f'{dir2}/{i}'
        out = f'{dir3}/{prefix}.afa'
        commands.append(f'mafft --quiet {_in} >{out}')

    print('running mafft...')
    with Pool(processes=threads) as p:
        with tqdm(total=len(commands)) as pbar:
            for i in p.imap(os.system, commands):
                pbar.update()

def build_tree(indir):
    'build phylogenetic tree by iqtree'
    outdir = f'{indir}/iqtree'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    command = ['iqtree', '-quiet']
    command.extend(['-s', f'{indir}/concatenated.afa'])
    command.extend(['-spp', f'{indir}/partition.txt'])
    command.extend(['-pre', f'{outdir}/iqtree'])
    command.extend(['-m', 'MFP'])
    print('building phylogenetic tree...')
    os.system(' '.join(command))

def main():
    args = parse_args()
    outdir=args.outdir
    script_path = os.path.abspath(__file__)
    data_dir = os.path.dirname(script_path) + '/data'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    dir1 = f'{outdir}/marker_gene'
    dir2 = f'{outdir}/gene_fasta'
    dir3 = f'{outdir}/aligned_fasta'
    dir4 = f'{outdir}/concatenated'
    tmp = f'{outdir}/tmp'
    for subdir in [dir1, dir2, dir3, dir4, tmp]:
        if not os.path.exists(subdir):
            os.mkdir(subdir)
    
    marker_genes = extract_marker_gene(args.genome_dir, dir1, tmp, data_dir, args.threads, args.force)
    a_gene_a_fasta(dir1, dir2)
    align(dir2, dir3, args.threads)
    concatenate(dir3, dir4)
    build_tree(dir4)

if __name__ == '__main__':
    main()
