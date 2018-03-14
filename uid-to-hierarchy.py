import ArgumentParser
import Entrez
import pandas as pd

def main():

    args = get_args()

    uid_df = pd.read_csv(args.uid, sep='\t', header=None, columns=['qid', 'uid', 'e'])
    uniq_uids = list(set(uid_df['uid']))
    uid_hier_df = get_hier(uniq_uids)

def get_hier(uids):
    '''
    Get taxonomic hierarchy of each UID
    '''



    return uid_hier_df

def get_args():
    '''
    Get command line arguments
    '''

    parser = ArgumentParser()
    parser.add_argument('-u', '--uid', help='Path to DIAMOND NCBI Taxonomy UID output file')
    parser.add_argument('-o', '--out', help='Path to taxonomic hierarchy table output')

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    main()