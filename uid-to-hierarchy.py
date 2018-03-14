import ArgumentParser
from Bio import Entrez
import pandas as pd
import sys
import time

from collections import OrderedDict
from functools import partial
from multiprocessing import current_process, Pool

Entrez.email = 'samuelmiller10@gmail.com'
target_ranks = set([
    'species', 
    'genus', 
    'family', 
    'order', 
    'class', 
    'phylum', 
    'superkingdom'
])

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

    uid_df = pd.read_csv(args.uid, sep='\t', header=None, columns=['qid', 'uid', 'e'])
    assert pd.api.types.is_numeric_dtype(uid_df['uid'])
    # e-values should include "e" for exponentiation
    assert pd.api.types.is_object_dtype(uid_df['e'])

    uniq_uids = list(set(uid_df['uid']))
    # 0 is the value when no alignment is found
    uniq_uids.remove(0)
    one_pct_tot = len(uniq_uids) / 100 / threads
    procedure = 'Taxonomic hierarchy recovery'
    rank_record = OrderedDict().fromkeys(target_ranks)
    mp_pool = Pool(threads)
    hiers = mp_pool.map(get_hier, uniq_uids)
    mp_pool.close()
    mp_pool.join()

def get_hier(uid, rank_record):
    '''
    Get taxonomic hierarchy from UID
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
                print('Waiting for Entrez to process UID ', str(uid), flush=True)
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
        for rank in target_ranks.difference(hier_ranks):
            rank_dict[rank] = ''

        hier = [rank for rank in rank_dict.values()]
        if hier[-1] not in superkingdoms:
            print(
                'Querying Entrez again with UID: ', str(uid), ' due to corrupted return data', 
                flush=True
            )
        else:
            break

    return hier

def get_lineage(uid):

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
    if config.verbose[0]:
        # Clear to end of line
        sys.stdout.write('\033[K')
        sys.stdout.write(msg + '\r')
        sys.stdout.flush()

def get_args():
    '''
    Get command line arguments
    '''

    parser = ArgumentParser()
    parser.add_argument('-u', '--uid', help='Path to DIAMOND NCBI Taxonomy UID output file')
    parser.add_argument('-o', '--out', help='Path to taxonomic hierarchy table output')
    parser.add_argument('-t', '--thread', help='Number of threads to use')

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    main()