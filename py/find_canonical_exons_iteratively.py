#!/usr/bin/env python3
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from statistics import median
from scipy.signal import find_peaks
from scipy.cluster import hierarchy

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to PAF file of read alignments")
    parser.add_argument("-i",
                        "--iterations",
                        type=int,
                        default=2,
                        help="Number of iterations for segmentation")
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
        if len(rid_to_intervals) > 99:
            break
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

def get_hetero_rids(band, low=0.1, high=0.9):
    return [rid for rid,cov in enumerate(band) if cov < high and cov > low]

def merge_peaks(peaks_a, peaks_b, range_len):
    result = list()
    peaks = sorted(set([p for p in peaks_a] + [p for p in peaks_b]))
    idx_s = 0
    while idx_s < len(peaks):
        for idx_e in range(idx_s+1,len(peaks)+1):
            if idx_e == len(peaks) or peaks[idx_e] - peaks[idx_e-1] > range_len:
                break
        result.append(int(median(peaks[idx_s:idx_e])))
        idx_s = idx_e
    return result

def find_exons(pos_to_rid, interval, rids, range_len=15):
    rids = set(rids)
    print('Finding exons in interval [{}:{}] with {} reads'.format(interval[0], interval[1], len(rids)))
    N = len(rids)
    x = list()
    y_rmved = list()
    y_added = list()

    for i in range(interval[0], interval[1]):
        j = max(interval[0], i-range_len)
        i_rids = pos_to_rid[i] & rids
        j_rids = pos_to_rid[j] & rids

        x.append(i)
        y_rmved.append(len(j_rids - i_rids))
        y_added.append(len(i_rids - j_rids))
    peaks_rmv, _ = find_peaks(y_rmved, height =max(N*0.05,5), distance=range_len, prominence=N*0.05)
    peaks_rmv    = [i+interval[0] for i in peaks_rmv]
    peaks_add, _ = find_peaks(y_added, height =max(N*0.05,5), distance=range_len, prominence=N*0.05)
    peaks_add    = [i+interval[0] for i in peaks_add]

    peaks = merge_peaks(peaks_a=peaks_rmv, peaks_b=peaks_add, range_len=range_len)
    exons = [interval]
    for idx,peak in enumerate(peaks):
        exons.append([peak,exons[-1][1]])
        exons[-2][1] = peak
    return exons

def plot_data(data, coverage, rid_to_intervals, N, M, out_path):
    L = len(data)
    cmap = plt.get_cmap('Greys')
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cmap_ins = plt.get_cmap('gnuplot2_r')
    norm_ins = mpl.colors.Normalize(vmin=10, vmax=800)
    top = 0.95
    bottom = 0.05
    step = (top-bottom)/N

    fig, axes = plt.subplots(L, 2, sharex='col', sharey='row', figsize=(50,8*L), squeeze=False)
    fig.suptitle('N = {} M = {} L = {}'.format(N,M,L))
    for axis_id,(exons,matrix) in enumerate(data):
        ax0 = axes[axis_id][0]
        ax1 = axes[axis_id][1]
        ax0.set_title('Stage {} with {} exons. Average coverage per exon'.format(axis_id, len(exons)))
        ax1.set_title('Stage {} with {} exons. True coverage'.format(axis_id, len(exons)))
        print('matrix.shape {} L = {} len(exons) {}'.format(matrix.shape, axis_id, len(exons)))
        reads_order = hierarchy.leaves_list(hierarchy.linkage(matrix, 'single'))
        for read_count,rid in enumerate(reads_order):
            if read_count % 50 == 0:
                print(read_count)
            h = top-step*read_count
            # Plot canonical exon boundies
            for start,end in exons[:-1]:
                ax0.plot([end,end], [0,1], color='blue', linestyle='dashed', zorder=100)
                ax1.plot([end,end], [0,1], color='blue', linestyle='dashed', zorder=100)
            # Plot average coverage of the read per exon
            for bid,(start,end) in enumerate(exons):
                print(start,end)
                ax0.plot([start,end], [h,h], color=cmap(norm(matrix[rid][bid])), lw=0.25, zorder=0)
            # Plot full alignments of the read
            for idx in range(len(rid_to_intervals[rid])):
                (cur_t_start,cur_t_end), (cur_q_start,cur_q_end) = rid_to_intervals[rid][idx]
                ax1.plot([cur_t_start,cur_t_end], [h,h], color='black', lw=0.25, zorder=0, alpha=0.5)
                if idx - 1 >= 0:
                    (_,_), (_,lst_q_end)   = rid_to_intervals[rid][idx-1]
                else:
                    lst_q_end   = 0
                dist_to_lst = cur_q_start - lst_q_end
                ax1.scatter(x=cur_t_start, y=h, s=0.25, color=cmap_ins(norm_ins(dist_to_lst)), zorder=5)
                if idx + 1 < len(rid_to_intervals[rid]):
                    (_,_), (nxt_q_start,_) = rid_to_intervals[rid][idx+1]
                else:
                    nxt_q_start = exons[-1][1]
                dist_to_nxt = nxt_q_start - cur_q_end
                ax1.scatter(x=cur_t_end,   y=h, s=0.25, color=cmap_ins(norm_ins(dist_to_nxt)), zorder=5)
        # Plot total coverage
        ax0.plot(range(len(coverage)), coverage, color='green', zorder=2.50)
        ax1.plot(range(len(coverage)), coverage, color='green', zorder=50)
    plt.savefig(out_path)


def main():
    args = parse_args()

    pos_to_rid,rid_to_intervals = read_paf(args.paf)
    N = len(rid_to_intervals)
    M = len(pos_to_rid)
    coverage = [len(rids)/N for rids in pos_to_rid]
    exons = [[0,M]]
    matrix = get_banded_matrix(N=N, pos_to_rid=pos_to_rid, intervals=exons)

    range_len = 15
    data = list()
    print('There are {} reads and {} positions'.format(N, M))
    for stage in range(args.iterations):
        old_exons = exons
        print('Running stage {} with {} exons and matrix shape {}'.format(stage, len(old_exons), matrix.shape))
        exons = list()
        for eid,interval in enumerate(old_exons):
            hetero_rids = get_hetero_rids(matrix[:,eid])
            print('Exon {} has {} hetero reads'.format(eid, len(hetero_rids)))
            if len(hetero_rids) > min(N*0.05,5) and interval[1]-interval[0] > range_len:
                exons.extend(find_exons(pos_to_rid=pos_to_rid, interval=interval, rids=hetero_rids))
            else:
                exons.append(interval)
        print('After re-segmentation there are {} exons'.format(len(exons)))
        matrix = get_banded_matrix(N=N, pos_to_rid=pos_to_rid, intervals=exons)
        data.append((exons, np.copy(matrix)))
    plot_data(data=data, coverage=coverage, rid_to_intervals=rid_to_intervals, N=N, M=M, out_path='{}.pdf'.format(args.out_prefix))

if __name__ == "__main__":
    main()
