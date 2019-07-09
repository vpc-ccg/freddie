#!/usr/bin/env python3
import argparse
import random
import numpy as np
import matplotlib as mpl
from scipy.cluster import hierarchy
mpl.use('Agg')
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    parser.add_argument("-t",
                        "--tsv",
                        type=str,
                        required=True,
                        help="Path to TSV file")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to PAF file of read alignments")
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

def get_banded_matrix(N, pos_to_rid, cut_points):
    height = N
    width = len(cut_points)-1
    banded_matrix = np.zeros(shape=(height, width), dtype=float)
    for rid in range(height):
        for bid in range(width):
            start = cut_points[bid]
            end = cut_points[bid+1]
            l = [rid in pos_to_rid[pos] for pos in range(start,end)]
            banded_matrix[rid][bid] = np.mean(l)
    return banded_matrix

def plot_reads(N, coverage, matrix, peaks, reads_order, rid_to_intervals, out_path):

    fig, axes = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(30,20))
    ax1 = axes[0]
    ax2 = axes[1]
    ax1.set_title('N = {}'.format(N))
    cmap = plt.get_cmap('Greys')
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cmap_ins = plt.get_cmap('gnuplot2_r')
    norm_ins = mpl.colors.Normalize(vmin=10, vmax=800)
    cax = plt.axes([0.95, 0.15, 0.01, 0.3])

    mpl.colorbar.ColorbarBase(norm=norm_ins, cmap=cmap_ins,orientation='vertical', ax=cax)
    top = 0.95
    bottom = 0.05
    step = (top-bottom)/N
    for idx,rid in enumerate(reads_order):
        if idx % 50 == 0:
            print(idx)
        h = top-step*idx
        for bid,val in enumerate(matrix[rid]):
            start = peaks[bid]
            end   = peaks[bid+1]
            ax1.plot([start,end], [h,h], color=cmap(norm(val)), lw=0.25, zorder=0)
        for idx in range(len(rid_to_intervals[rid])):
            (cur_t_start,cur_t_end), (cur_q_start,cur_q_end) = rid_to_intervals[rid][idx]
            ax2.plot([cur_t_start,cur_t_end], [h,h], color='black', lw=0.25, zorder=0, alpha=0.5)
            if idx - 1 >= 0:
                (_,_), (_,lst_q_end)   = rid_to_intervals[rid][idx-1]
            else:
                lst_q_end   = 0
            dist_to_lst = cur_q_start - lst_q_end
            ax2.scatter(x=cur_t_start, y=h, s=0.25, color=cmap_ins(norm_ins(dist_to_lst)), zorder=5)

            if idx + 1 < len(rid_to_intervals[rid]):
                (_,_), (nxt_q_start,_) = rid_to_intervals[rid][idx+1]
            else:
                nxt_q_start = peaks[-1]
            dist_to_nxt = nxt_q_start - cur_q_end
            ax2.scatter(x=cur_t_end,   y=h, s=0.25, color=cmap_ins(norm_ins(dist_to_nxt)), zorder=5)
    ax1.plot(range(len(coverage)), coverage, color='green', zorder=50)
    ax2.plot(range(len(coverage)), coverage, color='green', zorder=2.50)
    for peak in peaks:
        ax1.plot([peak,peak], [0,1], color='blue', linestyle='dashed', zorder=100)
        ax2.plot([peak,peak], [0,1], color='blue', linestyle='dashed', zorder=100)
    plt.savefig(out_path)

def main():
    args = parse_args()

    pos_to_rid,rid_to_intervals = read_paf(args.paf)
    N = len(rid_to_intervals)
    coverage = [len(rids)/N for rids in pos_to_rid]
    peaks = [int(p) for p in open(args.peaks).readlines()]
    banded_matrix = get_banded_matrix(N=N, pos_to_rid=pos_to_rid, cut_points=peaks)
    print(banded_matrix.shape)
    Z = hierarchy.linkage(banded_matrix, 'single')
    leaves = hierarchy.leaves_list(Z)
    out_file = open('{}.leaves.txt'.format(args.out_prefix), 'w+')
    for leaf in leaves:
        print(leaf, file=out_file)
    out_file.close()
    out_path = '{}.pdf'.format(args.out_prefix)
    plot_reads(N=N, coverage=coverage, matrix=banded_matrix, peaks=peaks, rid_to_intervals=rid_to_intervals, reads_order=leaves, out_path=out_path)

if __name__ == "__main__":
    main()
