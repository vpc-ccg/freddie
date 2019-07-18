#!/usr/bin/env python3
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean,stdev,median
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
    parser.add_argument("-op",
                        "--out-prefix",
                        type=str,
                        required=True,
                        help="Output prefix that does not include .EXT part")
    args = parser.parse_args()
    return args

def jaccard(a,b):
    i = len(a & b)
    j = len(a | b)
    if j == 0:
        return 0
    else:
        return i/j

def plot_coverage(coverage, transcripts, pos_to_rid, out_path, N):
    f, (ax2, ax3, ax1) = plt.subplots(3, 1, figsize=(30,5*(3+1)), sharex=True)
    ax1.plot(range(len(coverage)), coverage)
    ax3.plot(range(len(coverage)), coverage)
    M = len(transcripts)
    h_step = 0.8/M
    for idx,(exons,introns) in enumerate(transcripts):
        h = 0.9 - idx*h_step
        for exon in exons:
            ax1.plot(exon, [h,h], color='black', marker='o', alpha=0.8)
            ax2.plot(exon, [h*N,h*N], color='black', marker='o', alpha=0.8)
            ax3.plot(exon, [h,h], color='black', marker='o', alpha=0.8)
        for intron in introns:
            ax1.plot(intron, [h,h], color='gray', alpha=0.2)
            ax2.plot(intron, [h*N,h*N], color='gray', alpha=0.2)
            ax3.plot(intron, [h,h], color='gray', alpha=0.2)
    neighborhood_size = 15
    x_jaccard = list()
    y_jaccard_r = list()
    y_jaccard_l = list()
    plt.title('N = {}'.format(N))
    for i in range(len(pos_to_rid)):
        if len(pos_to_rid[i])/N < 0.01 or len(pos_to_rid[i]) < 3:
            x_jaccard.append(i)
            y_jaccard_r.append(0.0)
            y_jaccard_l.append(0.0)
            continue
        neighbors_r = set()
        prev = set()
        range_start = max(i-neighborhood_size, 0)
        range_end = i
        for j in range(range_start, range_end):
            prev = prev | pos_to_rid[j]
        y_jaccard_r.append(jaccard(prev, pos_to_rid[i]))

        neighbors_l = set()
        prev = set()
        range_start = i+1
        range_end = min(i+neighborhood_size+1, len(pos_to_rid))
        for j in range(range_start, range_end):
            prev = prev | pos_to_rid[j]
        y_jaccard_l.append(jaccard(prev, pos_to_rid[i]))
        x_jaccard.append(i)

    ax1.plot(x_jaccard, y_jaccard_r, color='r', alpha=0.25)
    ax1.plot(x_jaccard, y_jaccard_l, color='g', alpha=0.25)

    x_exonic_pos = list()
    y_exonic_jac_r = list()
    y_exonic_jac_l = list()

    cmap_r = plt.get_cmap('Reds')
    cmap_l = plt.get_cmap('Greens')
    for i in range(len(pos_to_rid)+1):
        if i == len(pos_to_rid) or len(pos_to_rid[i])/N < 0.01 or len(pos_to_rid[i]) < 3:
            if len(x_exonic_pos) >= 30:
                mu_r = mean(y_exonic_jac_r)
                sigma_r = stdev(y_exonic_jac_r)
                mu_l = mean(y_exonic_jac_l)
                sigma_l = stdev(y_exonic_jac_l)
                print(x_exonic_pos[0], x_exonic_pos[-1])
                stats_text = 'μ_r: {:.2f} σ_r: {:.2f} | μ_l: {:.2f} σ_l: {:.2f}'.format(mu_r, sigma_r, mu_l, sigma_l)
                print(stats_text)
                ax1.plot([x_exonic_pos[0],x_exonic_pos[-1]], [1.05,1.05], color='black', marker='x', alpha=1)
                text_x=x_exonic_pos[0]+(x_exonic_pos[-1]-x_exonic_pos[0])/10
                ax1.text(x=text_x, y=1.07, s=stats_text)

                norm_r = mpl.colors.Normalize(vmin=sigma_r, vmax=3*sigma_r)
                norm_l = mpl.colors.Normalize(vmin=sigma_l, vmax=3*sigma_l)

                for x,y_r,y_l in zip(x_exonic_pos,y_exonic_jac_r,y_exonic_jac_l):
                    if sigma_r > 0:
                        deviation_r = abs(mu_r-y_r)/sigma_r
                        if deviation_r > 0.5:
                            ax1.scatter(x=x,y=y_r, color=cmap_r(norm_r(deviation_r)), marker='>')
                    if sigma_l > 0:
                        deviation_l = abs(mu_l-y_l)/sigma_l
                        if deviation_l > 0.5:
                            ax1.scatter(x=x,y=y_l, color=cmap_l(norm_l(deviation_l)), marker='<')
                x_exonic_pos = list()
                y_exonic_jac_r = list()
                y_exonic_jac_l = list()
            continue
        x_exonic_pos.append(i)
        y_exonic_jac_r.append(y_jaccard_r[i])
        y_exonic_jac_l.append(y_jaccard_l[i])

    y_rmv_r = list()
    y_add_r = list()
    for i,rids in enumerate(pos_to_rid):
        j_r = max(0, i-neighborhood_size)
        rids_rmved_r = len(pos_to_rid[j_r] - rids)
        rids_added_r = len(rids - pos_to_rid[j_r])
        y_rmv_r.append(rids_rmved_r)
        y_add_r.append(rids_added_r)
    ax2.plot(range(len(pos_to_rid)), y_rmv_r, color='#e41a1c', alpha=0.5)
    ax2.plot(range(len(pos_to_rid)), y_add_r, color='#377eb8', alpha=0.5)

    peaks_rmv, _ = find_peaks(y_rmv_r, height=max(N*0.01,3), distance=neighborhood_size, prominence=N*0.025)
    peaks_add, _ = find_peaks(y_add_r, height=max(N*0.01,3), distance=neighborhood_size, prominence=N*0.025)
    ax2.plot(peaks_rmv, [y_rmv_r[i] for i in peaks_rmv], "x", color='#e41a1c')
    ax2.plot(peaks_add, [y_add_r[i] for i in peaks_add], "x", color='#377eb8')

    peaks_final = [0]
    peaks_dense = []
    peaks_all = sorted(set([p for p in peaks_rmv] + [p for p in peaks_add] + [0,len(pos_to_rid)-1]))
    idx = 1
    for idx in range(1, len(peaks_all)-1):
        prv_peak = peaks_all[idx-1]
        cur_peak = peaks_all[idx]
        nxt_peak = peaks_all[idx+1]
        before_dist = abs(cur_peak - prv_peak)
        after_dist  = abs(cur_peak - nxt_peak)
        if before_dist > neighborhood_size*2 and after_dist > neighborhood_size*2:
            peaks_final.append(cur_peak)
        else:
            peaks_dense.append(cur_peak)
    peaks_final.append(len(pos_to_rid)-1)

    peaks_dense_final = list()
    idx_s = 0
    for idx,peak in enumerate(peaks_dense):
        print(peak)
        idx_e = idx
        if idx+1 == len(peaks_dense) or peaks_dense[idx+1] - peak > neighborhood_size*2:
            print('d: ',peaks_dense[idx_s:idx_e+1])
            peaks_dense_final.append(int(median(peaks_dense[idx_s:idx_e+1])))
            idx_s = idx+1
    print('A:', peaks_all)
    print('F:', peaks_final)
    print('D:', peaks_dense)
    print('DF:', peaks_dense_final)
    for peak in peaks_final:
        ax3.plot([peak,peak], [0,1], color='black', linestyle='solid', lw=0.75)
    for peak in peaks_dense_final:
        ax3.plot([peak,peak], [0,1], color='black', linestyle='dashed', lw=0.75)
    for peak in peaks_dense:
        ax3.plot([peak,peak], [0,1], color='gray', linestyle='dotted', lw=0.75, alpha=0.50)

    plt.tight_layout()
    plt.savefig(out_path)
    return sorted(peaks_final+peaks_dense_final)

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
            rid_to_intervals[rid].append([t_start, t_end])
            for i in range(t_start, t_end):
                pos_to_rid[i].add(rid)
    for intervals in rid_to_intervals.values():
        intervals.sort()
    return pos_to_rid,rid_to_intervals

def get_tsv_ticks(tsv):
    transcripts = list()
    starts = list()
    ends = list()
    for line in open(tsv):
        exons = list()
        introns = list()
        line = line.rstrip().split('\t')
        for interval in line[3].split(','):
            interval = interval.split('-')
            start = int(interval[0])
            end = int(interval[1])
            if len(exons) > 0:
                introns.append((exons[-1][1], start))
            exons.append((start,end))
        transcripts.append((exons,introns))
    return transcripts

def main():
    args = parse_args()

    transcripts = get_tsv_ticks(args.tsv)
    pos_to_rid,rid_to_intervals = read_paf(args.paf)
    N = len(rid_to_intervals)
    coverage = [len(rids)/N for rids in pos_to_rid]
    outpath = '{}.pdf'.format(args.out_prefix)
    peaks = plot_coverage(coverage=coverage, transcripts=transcripts, pos_to_rid=pos_to_rid, out_path=outpath, N=N)
    outpath = '{}.txt'.format(args.out_prefix)
    out_file = open(outpath, 'w+')
    for i in peaks:
        print(i, file=out_file)
    out_file.close()

if __name__ == "__main__":
    main()
