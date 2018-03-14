import ArgumentParser
from Bio import Entrez
import pandas as pd
import sys

from functools import partial
from multiprocessing import current_process, Pool

Entrez.email = 'samuelmiller10@gmail.com'

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
    one_pct_tot = len(uniq_uids) / 100 / threads
    procedure = 'Taxonomic hierarchy recovery'
    mp_pool = Pool(threads)
    hiers = mp_pool.map(uniq_uids)
    mp_pool.close()
    mp_pool.join()

    uid_hier_df = get_hier(uniq_uids)

def get_hier(uid):
    '''
    Get taxonomic hierarchy from UID
    '''

    mp_pool = Pool(threads)


    return uid_hier_df

def get_lineage(uid):

def print_prog(procedure, one_pct_tot):
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