import argparse
from gurobipy import *

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    parser.add_argument("-d",
                        "--data-matrix",
                        type=str,
                        required=True,
                        help="Path to DATA file of reads by exons matrix")
    parser.add_argument("-op",
                        "--out_prefix",
                        type=str,
                        required=True,
                        help="Output prefix that does not include .EXT part")
    args = parser.parse_args()
    return args

def read_matrix(data_matrix):
    result = list()
    for line in open(data_matrix):
        line = line.rstrip()
        result.append([int(x) for x in line])
    return result

def main():
    args = parse_args()

    R = read_matrix(data_matrix=args.data_matrix)
    C = list()
    out_file = open('{}.txt'.format(args.out_prefix), 'w+')
    for exons in R:
        C.append([0 for _ in range(len(exons))])
        try:
            first_one = next((i for i, x in enumerate(exons) if x ==1))
            last_one  = len(exons)-1 - next((i for i, x in enumerate(reversed(exons)) if x ==1))
        except StopIteration:
            continue
        for i in range(first_one, last_one+1):
            C[-1][i] = 1
        out_file.write(''.join([str(x) for x in exons]))
        out_file.write('\n')
        out_file.write(''.join([str(x) for x in C[-1]]))
        out_file.write('\n')

    out_file.close()



if __name__ == "__main__":
    main()
