#!/usr/bin/env python3
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as pp
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to PAF file of read alignments")
    parser.add_argument("-t",
                        "--tsv",
                        type=str,
                        required=True,
                        help="Path to TSV file")
    parser.add_argument("-s",
                        "--sigma",
                        type=float,
                        default=5.0,
                        help="Sigma value for gaussian_filter1d")
    parser.add_argument("-op",
                        "--out-prefix",
                        type=str,
                        required=True,
                        help="Output prefix that does not include .TXT part")
    args = parser.parse_args()
    return args

def read_paf(paf, range_len=5):
    is_first = True
    pos_to_rid = list()
    read_name_to_id = dict()
    rid_to_intervals = dict()
    for line in open(paf):
        line = line.rstrip().split('\t')
        if is_first:
            t_len = int(line[6])
            t_name = line[5]
            is_first = False
            pos_to_rid = [set() for _ in range(t_len)]
        if t_len != int(line[6]) or t_name != line[5]:
            print("Multiple targets detected in PAF file!", file=stderr)
            print(line, file=stderr)
            exit(-1)
        name = line[0]
        if not name in read_name_to_id:
            rid = len(read_name_to_id)
            read_name_to_id[name] = rid
            rid_to_intervals[rid] = list()
        rid = read_name_to_id[name]
        if any('oc:c:1' in tag for tag in line[12:]):
            t_start = max(0, int(line[7]) - 1)
            t_end = int(line[8]) + 1
            q_start = max(0, int(line[2]) - 1)
            q_end = int(line[3]) + 1
            t_interval = (t_start, t_end)
            q_interval = (q_start, q_end)
            rid_to_intervals[rid].append((t_interval, q_interval))
    for intervals in rid_to_intervals.values():
        intervals.sort()
    for rid,intervals in rid_to_intervals.items():
        new_intervals = list()
        for idx,(t_interval,q_interval) in enumerate(intervals):
            if idx == 0:
                new_intervals.append((t_interval,q_interval))
                continue
            (t_interval_prv,q_interval_prv) = new_intervals[-1]
            if t_interval[0] - t_interval_prv[1] < range_len and q_interval_prv[0] - q_interval_prv[1] < range_len:
                new_intervals[-1] = (
                    (t_interval_prv[0],t_interval[1]),
                    (q_interval_prv[0],q_interval[1]),
                )
            else:
                new_intervals.append((t_interval,q_interval))
        rid_to_intervals[rid] = new_intervals
    for rid,intervals in rid_to_intervals.items():
        for (t_start, t_end),(_, _) in intervals:
            for i in range(t_start, t_end):
                pos_to_rid[i].add(rid)
    return pos_to_rid,rid_to_intervals,read_name_to_id

def get_brks(tsv):
    brks = set()
    for line in open(tsv):
        line = line.rstrip().split()
        for i in line[3].split(','):
            brks.add(int(i.split('-')[0]))
            brks.add(int(i.split('-')[1]))
    return sorted(brks)

def get_splice(gene_len, rid_to_intervals):
    X = [x for x in range(0,gene_len+1)]
    Y_i = np.zeros(len(X))
    Y_o = np.zeros(len(X))
    for rid,intervals in rid_to_intervals.items():
        for (t_start, t_end),(_, _) in intervals:
            Y_i[t_start] =+1
            Y_o[t_end]   =+1
    Y_a = Y_i + Y_o
    return X, Y_i, Y_o, Y_a

def plot(Y_roll, peaks, brks, outpath):
    pp.figure(figsize=(20,5))
    pp.plot(Y_roll)
    pp.plot(peaks, [Y_roll[p] for p in peaks], "x")
    for b in brks:
        pp.axvline(b,ymin=0, ymax=max(Y_roll)*1.05, color='green', ls='dashed', alpha=0.4, linewidth=1)
    pp.title('|p|={}'.format(len(peaks)))
    pp.savefig(fname=outpath)

def main():
    args = parse_args()

    pos_to_rid,rid_to_intervals,read_name_to_id = read_paf(args.paf)
    brks = get_brks(args.tsv)
    if len(rid_to_intervals)>0:
        gene_len = int(open(args.paf).readline().split('\t')[6])
    else:
        gene_len = 1
    X, Y_i, Y_o, Y_a = get_splice(gene_len=gene_len, rid_to_intervals=rid_to_intervals)

    if not args.sigma > 0.0:
        print('Sigma is not positive float. Not doing any filtering')
        Y_roll=Y_a
    else:
        Y_roll=gaussian_filter1d(Y_a,args.sigma)
    peaks, _ = find_peaks(Y_roll)
    plot(Y_roll=Y_roll, peaks=peaks, brks=brks, outpath='{}.pdf'.format(args.out_prefix))
    out_file = open('{}.txt'.format(args.out_prefix), 'w+')
    for i in peaks:
        print(i, file=out_file)
    out_file.close()

if __name__ == "__main__":
    main()
