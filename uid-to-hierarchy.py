import ArgumentParser

def main():

    args = get_args()

def get_args():
    '''
    Get command line arguments
    '''

    parser = ArgumentParser()
    parser.add_argument('-u', '--uids', help='Path to DIAMOND NCBI Taxonomy id output file')
    parser.add_argument('-o', '--out', help='Path to taxonomic hierarchy table output')

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    main()