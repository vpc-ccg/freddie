#!/usr/bin/env python3
import argparse
from sys import stderr
from sklearn.cluster import DBSCAN
from sklearn.metrics.cluster import adjusted_rand_score
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster all barcodes with hamming distance thershold")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to PAF file")
    parser.add_argument("-m",
                        "--min-samples",
                        type=str,
                        default=2,
                        help="sklearn DBSCAN min_samples param")
    parser.add_argument("-e",
                        "--eps",
                        type=str,
                        default=1,
                        help="sklearn DBSCAN eps param")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        help="Put TSV file with results of clustering")
    args = parser.parse_args()
    return args

def read_paf(paf):
    reads = dict()
    is_first = True
    for line in open(paf):
        line = line.rstrip().split('\t')
        if is_first:
            t_len = int(line[6])
            t_name = line[5]
            is_first = False
        if t_len != int(line[6]) or t_name != line[5]:
            print("Multiple targets detected in PAF file!", file=stderr)
            print(line, file=stderr)
            exit(-1)
        name = line[0]
        if not name in reads:
            reads[name] = [0 for _ in range(t_len)]
        if any('oc:c:1' in tag for tag in line[12:]):
            t_start = int(line[7])
            t_end = int(line[8])
            for i in range(t_start, t_end+1):
                reads[name][i] = 1
    names = sorted(list(reads.keys()))
    matrix = np.array([
        reads[name] for name in names
    ])
    return names,matrix

def get_ari_score(names, labels, true_labels):
    labels = [-(idx+1) if label < 0 else label for idx,label in enumerate(labels)]
    return adjusted_rand_score(true_labels, labels)

def main():
    args = parse_args()
    names,matrix = read_paf(paf=args.paf)

    true_labels = {x.split('_')[0] for x in names}
    true_labels = {label:idx for idx,label in enumerate(true_labels)}
    true_labels = [true_labels[x.split('_')[0]] for x in names]

    out_file = open(args.output, 'w+')
    for eps in np.arange(1, 10000, 100):
        for min_samples in range(2,3):
            clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(matrix)
            ari_score = get_ari_score(names=names, true_labels=true_labels, labels=list(clustering.labels_))
            print('{}\t{}\t{}'.format(min_samples, eps, ari_score), file=out_file)
    out_file.close()

if __name__ == "__main__":
    main()
