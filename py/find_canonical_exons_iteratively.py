#!/usr/bin/env python3
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from statistics import median
from scipy.signal import find_peaks

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to PAF file of read alignments")
    parser.add_argument("-l",
                        "--leaves",
                        type=str,
                        required=True,
                        help="Path to leaves order TXT file")
    parser.add_argument("-pk",
                        "--peaks",
                        type=str,
                        required=True,
                        help="Path to peaks TXT file")
    parser.add_argument("-op",
                        "--out_prefix",
                        type=str,
                        required=True,
                        help="Output prefix that does not include .EXT part")
    args = parser.parse_args()
    return args

def read_paf(paf):
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
            for i in range(t_start, t_end):
                pos_to_rid[i].add(rid)
    for intervals in rid_to_intervals.values():
        intervals.sort()
    return pos_to_rid,rid_to_intervals

def get_banded_matrix(N, pos_to_rid, intervals):
    height = N
    width = len(intervals)
    banded_matrix = np.zeros(shape=(height, width), dtype=float)
    for rid in range(height):
        for bid in range(width):
            start,end = intervals[bid]
            l = [rid in pos_to_rid[pos] for pos in range(start,end)]
            banded_matrix[rid][bid] = np.mean(l)
    return banded_matrix

def main():
    args = parse_args()

    pos_to_rid,rid_to_intervals = read_paf(args.paf)
    N = len(rid_to_intervals)
    M = len(pos_to_rid)
    coverage = [len(rids)/N for rids in pos_to_rid]
    exons = [(0,M)]

    for round_idx in range(5):
        matrix = get_banded_matrix(N=N, pos_to_rid=pos_to_rid, intervals=exons)
        old_exons = exons
        exons = list()
        for bid,interval in enumerate(exons):
            if is_hetero(matrix[bid]):
                peaks = find_exons(pos_to_rid=pos_to_rid, interval=interval)

            else:
                exons.append(old_exons[bid])




if __name__ == "__main__":
    main()
