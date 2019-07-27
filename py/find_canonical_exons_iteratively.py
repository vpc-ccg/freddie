#!/usr/bin/env python3
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from statistics import median
from scipy.signal import find_peaks
from scipy.cluster import hierarchy
from os import remove

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
                        default=10,
                        help="Number of iterations for segmentation")
    parser.add_argument("-op",
                        "--out-prefix",
                        type=str,
                        required=True,
                        help="Output prefix that does not include .EXT part")
    args = parser.parse_args()
    return args

def read_paf(paf, range_len=15):
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
    peaks_rmv, _ = find_peaks(y_rmved, height =max(N*0.025,10), distance=range_len)
    peaks_rmv    = [max(0,i+interval[0]-range_len//2) for i in peaks_rmv]
    peaks_add, _ = find_peaks(y_added, height =max(N*0.025,10), distance=range_len)
    peaks_add    = [max(0,i+interval[0]-range_len//2) for i in peaks_add]

    peaks = merge_peaks(peaks_a=peaks_rmv, peaks_b=peaks_add, range_len=range_len)
    exons = [interval]
    for idx,peak in enumerate(peaks):
        exons.append([peak,exons[-1][1]])
        exons[-2][1] = peak
    return x, y_rmved, y_added, exons

def simpify_singal(matrix, low=0.1, high=0.9):
    result = ['' for _ in range(matrix.shape[0])]
    for rid in range(matrix.shape[0]):
        for pos in range(matrix.shape[1]):
            x = matrix[rid,pos]
            if   x >= high:
                result[rid] += '1'
            elif x <= low:
                result[rid] += '0'
            else:
                result[rid] += '2'
    return result

def binarize_matrix(matrix, cutoff=0.5):
    return [[1 if matrix[i,j]>cutoff else 0 for j in range(matrix.shape[1])] for i in range(matrix.shape[0])]

def plot_data(data, coverage, rid_to_intervals, final_exon_intervals, final_matrix, N, M, range_len, out_prefix):
    L = len(data)
    fig, axes = plt.subplots(L+1, 1, sharex='col', sharey='row', figsize=(30,8*(L+1)), squeeze=False)
    fig.suptitle('N = {} M = {} L = {}'.format(N,M,L))
    for stage_id,exons in enumerate(data):
        ax0 = axes[stage_id][0]
        for eid,exon in enumerate(exons):
            h = [0.3*N, 0.6*N, 0.9*N][eid%3]
            if exon['fixed']:
                ax0.vlines(exon['interval'], ymin=0, ymax=N, color='gray', linestyle='solid', lw=1, alpha=0.5)
            else:
                ax0.vlines(exon['interval'], ymin=0, ymax=N, color='green', linestyle='dashed', lw=1, alpha=0.5)
                ax0.plot(exon['sig_x'], exon['sig_ry'], color='#e41a1c', alpha=0.5)
                ax0.plot(exon['sig_x'], exon['sig_ay'], color='#377eb8', alpha=0.5)
            ax0.text(x=exon['interval'][0]+10, y=h, s='{}({})'.format(exon['interval'],exon['h_cnt']), fontsize=12, rotation=90)
    # simplified_signal = simpify_singal(final_matrix)
    reads_order = hierarchy.leaves_list(hierarchy.linkage(final_matrix, 'single'))

    cmap_ins = plt.get_cmap('gnuplot2_r')
    norm_ins = mpl.colors.Normalize(vmin=10, vmax=800)
    top    = N*0.95
    bottom = N*0.05
    step = (top-bottom)/N
    ax0 = axes[L][0]
    for read_count,rid in enumerate(reads_order):
        if read_count % 50 == 0:
            print('{}'.format(read_count))
        h = top-step*read_count
        # Plot full alignments of the read
        for idx in range(len(rid_to_intervals[rid])):
            (cur_t_start,cur_t_end), (cur_q_start,cur_q_end) = rid_to_intervals[rid][idx]
            ax0.plot([cur_t_start,cur_t_end], [h,h], color='black', lw=0.25, zorder=0, alpha=0.5)
            if idx - 1 >= 0:
                (_,_), (_,lst_q_end)   = rid_to_intervals[rid][idx-1]
            else:
                lst_q_end   = 0
            dist_to_lst = cur_q_start - lst_q_end
            ax0.scatter(x=cur_t_start, y=h, s=0.25, color=cmap_ins(norm_ins(dist_to_lst)), zorder=5)
            if idx + 1 < len(rid_to_intervals[rid]):
                (_,_), (nxt_q_start,_) = rid_to_intervals[rid][idx+1]
            else:
                nxt_q_start = data[0][0]['interval'][1]
            dist_to_nxt = nxt_q_start - cur_q_end
            ax0.scatter(x=cur_t_end,   y=h, s=0.25, color=cmap_ins(norm_ins(dist_to_nxt)), zorder=5)
        # print(simplified_signal[rid])
        # if read_count > 0 and simplified_signal[reads_order[read_count-1]] != simplified_signal[rid]:
        #     print(simplified_signal[rid])
        #     ax0.hlines([h], xmin=0, xmax=data[0][0]['interval'][1], color='orange', linestyle='dashed')
    for interval in final_exon_intervals:
        ax0.vlines(interval, ymin=0, ymax=N, color='black', linestyle='solid', lw=2, )
    ax0.plot(range(len(coverage)), [N*c for c in coverage], color='green', zorder=50)
    plt.savefig('{}.pdf'.format(out_prefix))

def output_last_matrix(binary_matrix, outpath,):
    out_file = open(outpath, 'w+')
    for i in range(len(binary_matrix)):
        for j in range(len(binary_matrix[i])):
            print(binary_matrix[i][j], end='', file=out_file)
        print(file=out_file)
    out_file.close()

def merge_exons(exons, range_len):
    final_exon_intervals = list()
    eid = 0
    while eid < len(exons):
        s,e = exons[eid]['interval']
        while e-s < range_len and eid+1 < len(exons):
            eid += 1
            _,e = exons[eid]['interval']
        final_exon_intervals.append([s,e])
        eid += 1
    return final_exon_intervals

def output_final_exons(final_exon_intervals, outpath):
    out_file = open(outpath, 'w+')
    for s,e in final_exon_intervals:
        print('{}\t{}'.format(s,e), file=out_file)
    out_file.close()

def get_stretches_of_zeros(l):
    result = list()
    s = len(l)
    for i,v in enumerate(l):
        if l[i] == 1:
            if s < i-1:
                result.append([s,i])
            s = i
    return result

def get_unaligned(exons, zeros, r_intervals):
    result = list()
    for left_eid,right_eid in zeros:
        left_exon, right_exon = exons[left_eid], exons[right_eid]
        for idx,((t_start,_), (_,q_end)) in enumerate(r_intervals):
            if t_start > left_exon[1]:
                break
            left_end = q_end
        for (_,t_end), (q_start,_) in r_intervals[idx:]:
            right_start = q_start
            if t_end > right_exon[0]:
                break
        read_dist = right_start - left_end
        result.append(read_dist)
    return result

def output_read_unaligned_values(exons, binary_matrix, rid_to_intervals, outpath):
    out_file = open(outpath, 'w+')
    for rid,r_intervals in rid_to_intervals.items():
        zeros = get_stretches_of_zeros(l=binary_matrix[rid])
        unaligned = get_unaligned(exons=exons, zeros=zeros, r_intervals=r_intervals)
        for (i,j),l in zip(zeros, unaligned):
            print('{}-{}-{}'.format(i,j,abs(l)), end='\t', file=out_file)
        print(file=out_file)
    out_file.close()

def main():
    args = parse_args()

    pos_to_rid,rid_to_intervals = read_paf(args.paf)
    N = len(rid_to_intervals)
    M = len(pos_to_rid)
    coverage = [len(rids)/N for rids in pos_to_rid]

    exons = [
        dict(
            interval=[0,M],
            fixed=False,
            h_cnt=0,
            sig_x=list(),
            sig_ry=list(),
            sig_ay=list(),
        )
    ]
    data = [exons]
    range_len = 15
    for d in data:
        print(d)
    for stage_id in range(args.iterations):
        old_exons = data[stage_id]
        new_exons = list()
        matrix = get_banded_matrix(N=N, pos_to_rid=pos_to_rid, intervals=[exon['interval'] for exon in old_exons])
        print('Running stage {} with {} exons and matrix shape {}'.format(stage_id, len(old_exons), matrix.shape))
        for eid,exon in enumerate(old_exons):
            if exon['interval'][1]-exon['interval'][0] > 250:
                hetero_rids = [rid for rid in range(len(rid_to_intervals))]
            else:
                hetero_rids = get_hetero_rids(matrix[:,eid])
            exon['h_cnt'] = len(hetero_rids)
            exon['fixed'] = len(hetero_rids) <= min(N*0.10,25) or exon['interval'][1]-exon['interval'][0] < range_len
            if not exon['fixed']:
                print(exon['interval'])
                x, y_rmved, y_added, new_exon_intervals = find_exons(pos_to_rid=pos_to_rid, interval=exon['interval'].copy(), rids=hetero_rids)
                print(exon['interval'])
                exon['sig_x']  = x
                exon['sig_ry'] = y_rmved
                exon['sig_ay'] = y_added
                for new_exon_interval in new_exon_intervals:
                    new_exons.append(
                        dict(
                            interval=new_exon_interval,
                            fixed=False,
                            h_cnt=0,
                            sig_x=list(),
                            sig_ry=list(),
                            sig_ay=list(),
                        )
                    )
            else:
                print('Exon {} {} is fixed'.format(eid, exon['interval']))
                new_exons.append(exon)
        print('After re-segmentation there are {} exons'.format(len(new_exons)))
        if len(new_exons) == len(old_exons):
            print('The number of exons has not changed. Breaking...')
            break
        data.append(new_exons)
    final_exon_intervals = merge_exons(exons=data[-1], range_len=range_len)
    output_final_exons(final_exon_intervals=final_exon_intervals, outpath='{}.tsv'.format(args.out_prefix))
    final_matrix = get_banded_matrix(N=N, pos_to_rid=pos_to_rid, intervals=final_exon_intervals)
    binary_matrix = binarize_matrix(final_matrix)
    output_last_matrix(binary_matrix=binary_matrix, outpath='{}.data'.format(args.out_prefix))
    output_read_unaligned_values(exons=final_exon_intervals, binary_matrix=binary_matrix, rid_to_intervals=rid_to_intervals, outpath='{}.zeros_unaligned.tsv'.format(args.out_prefix))
    # plot_data(data=data, coverage=coverage, rid_to_intervals=rid_to_intervals, final_matrix=final_matrix, final_exon_intervals=final_exon_intervals, N=N, M=M, range_len=range_len, out_prefix=args.out_prefix)



if __name__ == "__main__":
    main()
