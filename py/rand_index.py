import sys
import argparse

from sklearn.metrics.cluster import adjusted_rand_score

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculates Rand Index for a clusters file")
    parser.add_argument("-r",
                        "--reads",
                        type=str,
                        required=True,
                        help="Input read FASTA file")
    parser.add_argument("-i",
                        "--isoforms",
                        type=str,
                        required=True,
                        help="Input predicted isoforms file")
    parser.add_argument("-o",
                        "--output-accuracy-results",
                        type=str,
                        required=False,
                        help="Output file where ARI score is outputted. Default: stdout")
    args = parser.parse_args()
    return args

def choose_2(n):
    return n*(n-1)//2

def main():
    args = parse_args()

    rid_to_tcid = dict()
    rid_to_pcid = dict()

    rnom_to_id = dict()
    tcnom_to_id = dict()
    pcnom_to_id = dict()

    for line in open(args.reads, 'r'):
        if line[0] != '>':
            continue

        rnom = line.rstrip().split()[0][1:]
        rid = rnom_to_id.get(rnom, len(rnom_to_id))
        rnom_to_id[rnom] = rid

        tcnom = rnom.split('_')[0]
        tcid = tcnom_to_id.get(tcnom, len(tcnom_to_id))
        tcnom_to_id[tcnom] = tcid

        rid_to_tcid[rid] = tcid

    if (args.output_accuracy_results):
        results = open(args.output_accuracy_results, 'w+')
    else:
        results = sys.stdout

    for line in open(args.isoforms, 'r'):
        if line[0] == '#':
            continue
        line = line.rstrip().split('\t')

        pcnom = line[0]
        pcid = pcnom_to_id.get(pcnom, len(pcnom_to_id))
        pcnom_to_id[pcnom] = pcid

        rnom = line[1]
        rid = rnom_to_id[rnom]
        rid_to_pcid[rid] = pcid

    Y = []
    X = []
    lonely_read_cluster_counter = -1
    for rid in rid_to_tcid:
        Y.append(rid_to_tcid[rid])
        if rid in rid_to_pcid:
            X.append(rid_to_pcid[rid])
        else:
            X.append(lonely_read_cluster_counter)
            lonely_read_cluster_counter -= 1
    print(adjusted_rand_score(X, Y), file=results)

if __name__ == "__main__":
    main()
