'''
Create taxonomic tables for Anvio
'''

import argparse
from Bio import Entrez
import pandas as pd
import sys
import time

from collections import OrderedDict
from functools import partial
from multiprocessing import current_process, Pool

Entrez.email = 'samuelmiller10@gmail.com'
target_ranks = [
    'species', 
    'genus', 
    'family', 
    'order', 
    'class', 
    'phylum', 
    'superkingdom'
]
superkingdoms=[
    'Archaea', 
    'Bacteria',
    'Eukaryota',
    'Viruses'    
]
tax_out_cols = [
    'gene_callers_id', 
    't_phylum', 
    't_class', 
    't_order', 
    't_family', 
    't_genus', 
    't_species'
]

# Globals set in parent process
threads = 1
one_pct_tot = 0
procedure = ''
# Global exclusively set and accessed in child process
prog_count = 0

def main():

    global threads, one_pct_tot

    args = get_args()
    threads = args.threads

    uid_df = pd.read_csv(args.uid, sep='\t', header=None, names=['qid', 'uid', 'e'])
    assert pd.api.types.is_numeric_dtype(uid_df['uid'])
    assert pd.api.types.is_numeric_dtype(uid_df['e'])
    
    uids = uid_df['uid'].tolist()
    uniq_uids = list(set(uids))
    # 0 is the value when no alignment is found
    uniq_uids.remove(0)
    # "Cellular organisms"
    uniq_uids.remove(131567)
    one_pct_tot = len(uniq_uids) / 100 / threads
    procedure = 'Taxonomic hierarchy recovery'
    rank_record = OrderedDict().fromkeys(target_ranks)
    mp_pool = Pool(threads)
    froz_get_hier = partial(get_hier, rank_record=rank_record)
    uniq_hiers = mp_pool.map(froz_get_hier, uniq_uids)
    mp_pool.close()
    mp_pool.join()

    uid_hier_dict = OrderedDict().fromkeys(uniq_uids)
    uid_hier_dict[0] = ['' for rank in target_ranks]
    uid_hier_dict[131567] = ['' for rank in target_ranks]
    for i, uniq_uid in enumerate(uniq_uids):
        uid_hier_dict[uniq_uid] = uniq_hiers[i]
    hiers = []
    for i, uid in enumerate(uids):
        hiers.append(uid_hier_dict[uid])
    hier_df = pd.DataFrame(hiers, columns=target_ranks)
    out_df = pd.concat([uid_df[['qid']], hier_df[target_ranks[::-1]]], axis=1)

    if bool(args.gene_table):
        make_misc_tbl(out_df, args.gene_table, args.split_fasta)

    # Write taxonomy table
    tax_out_df = out_df.drop('superkingdom', axis=1)
    tax_out_df.columns = tax_out_cols
    tax_out_df.to_csv(args.out, sep='\t', index=False)

    return

def make_misc_tbl(gene_tax_tbl, gene_contig_tbl, split_fasta, out):
    '''
    Make table with which each rank can be imported into Anvio as a layer (via anvi-import-misc-data)
    '''

    # Record the number of splits per contig
    split_count_dict = OrderedDict()
    with open(split_fasta) as handle:
        for line in handle.readlines():
            if line[0] == '>':
            id_parts = line.split('_split_')
            contig_id = id_parts[0].lstrip('>')
            try:
                split_count_dict[contig_id] += 1
            except KeyError:
                split_count_dict[contig_id] = 1

    # Annotate the taxonomy of each contig
    gene_contig_df = pd.read_csv(
        gene_contig_tbl, sep='\t', header=0, usecols=['gene_callers_id', 'contig']
    )
    gene_contig_df = gene_contig_df.merge(gene_tax_tbl, how=left, on='gene_callers_id')
    gene_contig_gb = gene_contig_df.groupby('contig', as_index=False)
    contig_tax_dict = OrderedDict([(k, []) for k in ['contig', 'split_count'] + target_ranks])
    for contig, group in gene_contig_gb:
        contig_tax_dict['contig'] = contig
        contig_tax_dict['split_count'] = split_count_dict[contig]
        for rank in target_ranks:
            freq_tax = group[rank].value_counts().index.tolist()
            if len(freq_tax) == 1:
                contig_tax_dict[rank] = freq_tax[0]
            else:
                for tax in freq_tax:
                    if freq_tax != '':
                        contig_tax_dict[rank] = tax
    contig_tax_df = pd.DataFrame(contig_tax_dict)
    
    # Annotate the taxonomy of each split
    splits = []
    contigs = []
    for contig, split_count in split_count_dict.items():
        for i in range(1, split_count + 1):
            split_num = str(i)
            split = contig + '_split_' + (5 - len(split_num)) * '0' + str(split_num)
            splits.append(split)
            contigs.append(contig)
    split_tax_df = pd.DataFrame({'split': splits, 'contig': contigs})
    split_tax_df = split_tax_df.merge(contig_tax_df, how=left, on='contig')

    split_tax_df.drop('contig', axis=1, inplace=True)
    new_basename = os.path.splitext(os.path.basename(out))[0] + '.splits'
    ext = os.path.splitext(os.path.basename(out))[1]
    path = os.path.join(
        os.path.dirname(out), 
        new_basename + ext
    )
    split_tax_df.to_csv(path, sep='\t', index=False)

    return

def get_hier(uid, rank_record):
    '''
    Get taxonomic hierarchy for UID from Entrez
    '''

    print_prog()

    # Entrez occasionally returns corrupted data
    bad_rtn = True
    for retrials in range(5):
        while True:
            try:
                data = Entrez.read(Entrez.efetch(db='Taxonomy', id=uid))
                break
            except:
                print('Waiting for Entrez to process UID', str(uid), flush=True)
                time.sleep(2)
        uid_taxon = data[0]['ScientificName']
        uid_rank = data[0]['Rank']
        hier_data = data[0]['LineageEx']

        hier_ranks = set()
        if uid_rank in target_ranks:
            rank_record[uid_rank] = uid_taxon
            hier_ranks.add(uid_rank)
        # Record higher ranks
        for entry in hier_data:
            higher_rank = entry['Rank']
            if higher_rank in target_ranks:
                rank_record[higher_rank] = entry['ScientificName']
                hier_ranks.add(higher_rank)
        # Record unresolvable ranks
        for rank in set(target_ranks).difference(hier_ranks):
            rank_record[rank] = ''

        hier = [rank for rank in rank_record.values()]
        if hier[-1] not in superkingdoms:
            print(hier, flush=True)
            print(
                'Querying Entrez again with UID: ', str(uid), ' due to corrupted return data', 
                flush=True
            )
        else:
            break

    return hier

def print_prog():
    '''
    Print percent progress
    '''

    if current_process()._identity[0] % threads == 1:
        global prog_count
        prog_count += 1
        if int(prog_count % one_pct_tot) == 0:
            pct = int(prog_count / one_pct_tot)
            if pct <= 100:
                print_over_line('{0}: {1}%'.format(procedure, pct))
                print_over_line(procedure + ':' + str(pct) + '%')    

def print_over_line(msg):
    # Clear to end of line
    sys.stdout.write('\033[K')
    sys.stdout.write(msg + '\r')
    sys.stdout.flush()

def get_args():
    '''
    Get command line arguments
    '''

    parser = argparse.ArgumentParser(
        help=(
            'Generate taxonomy table for addition to contigs table '
            'and optionally taxonomy table for addition to profile as miscellaneous data'
        )
    )
    parser.add_argument('uid', help='Path to DIAMOND NCBI Taxonomy UID output file')
    parser.add_argument('out', help='Path to taxonomic hierarchy table output')
    parser.add_argument('--gene_table', help='Path to anvi-output-gene-calls output')
    parser.add_argument('--split_fasta', help='Path to anvi-export-contigs --splits-mode output')
    parser.add_argument(
        '--threads', 
        default=1, 
        type=int, 
        help='Number of threads to use'
    )

    args = parser.parse_args()
    if bool(args.gene_table) != bool(args.split_fasta):
        parser.error('Both gene_table and split_fasta must be specified')

    return args

if __name__ == "__main__":
    main()