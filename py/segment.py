#!/usr/bin/env python3
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as pp
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from multiprocessing import Pool

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
    parser.add_argument("-c",
                        "--threads",
                        type=int,
                        default=4,
                        help="Number of threads for multiprocessing")
    parser.add_argument("-s",
                        "--sigma",
                        type=float,
                        default=5.0,
                        help="Sigma value for gaussian_filter1d")
    parser.add_argument("-lt",
                        "--low-threshold",
                        type=float,
                        default=0.10,
                        help="Low threshold under which the read will be considered as not covering a segment. Default: 0.1")
    parser.add_argument("-ht",
                        "--high-threshold",
                        type=float,
                        default=0.90,
                        help="High threshold above which the read will be considered as covering a segment. Default: 0.9")
    parser.add_argument("-op",
                        "--out-prefix",
                        type=str,
                        required=True,
                        help="Output prefix that does not include .TXT part")
    parser.add_argument("-v",
                        "--variance-factor",
                        type=float,
                        default=3.0,
                        help="The sigma factor to fix a candidate peak. The threshold is set as > mean+3*variance_factor. Default 3.0")
    parser.add_argument("-hp",
                        "--max-peak-cnt-per-segment",
                        type=int,
                        default=50,
                        help="Maximum number of peaks allowed per segment")
    parser.add_argument("-lo",
                        "--min-read-support-outside",
                        type=int,
                        default=3,
                        help="Minimum reads support for splice site to support a breakpoint")
    args = parser.parse_args()
    return args

def read_paf(paf, range_len=0):
    is_first = True
    pos_to_rid = list()
    read_name_to_id = dict()
    rid_to_intervals = dict()
    rid_to_len = dict()
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
            rid_to_len[rid] = int(line[1])
        rid = read_name_to_id[name]
        if any('oc:c:1' in tag for tag in line[12:]):
            t_start = int(line[7])
            t_end   = int(line[8])
            q_start = int(line[2])
            q_end   = int(line[3])
            t_interval = (t_start, t_end)
            q_interval = (q_start, q_end)
            rid_to_intervals[rid].append((t_interval, q_interval))
            assert rid_to_len[rid] == int(line[1])
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
    return pos_to_rid,rid_to_len,rid_to_intervals,read_name_to_id,t_len

def annotated_breakpoints(tsv):
    """
    Get splicing breakpoints from annotation TSV file for plotting purposes only.
    """
    annotated_brks = set()
    for line in open(tsv):
        line = line.rstrip().split()
        for i in line[3].split(','):
            annotated_brks.add(int(i.split('-')[0]))
            annotated_brks.add(int(i.split('-')[1]))
    return sorted(annotated_brks)

def get_splice(gene_len, rid_to_intervals):
    """
    Get splicing in and out for each postion on the gene.
    """
    X = [x for x in range(0,gene_len+1)]
    Y_i = np.zeros(len(X))
    Y_o = np.zeros(len(X))
    for rid,intervals in rid_to_intervals.items():
        for (t_start, t_end),(_, _) in intervals:
            Y_i[t_start] =+1
            Y_o[t_end]   =+1
    Y_a = Y_i + Y_o
    return X, Y_i, Y_o, Y_a

def get_high_var_peaks(peaks, Y_roll, variance_factor):
    """
    Returns a set of peak indices that are above mean+variance_factor*std
    of the signal
    """
    thresh = Y_roll[Y_roll>0].mean() + variance_factor*Y_roll[Y_roll>0].std()
    peaks_vrnc = set()
    for idx,p in enumerate(peaks):
        if Y_roll[p] > thresh:
            peaks_vrnc.add(idx)
    return peaks_vrnc

def get_desert_bound_peaks(peaks, pos_to_rid, width=50):
    """
    Returns a list of peak indices of peaks that neighbor a region of deserts
    of size >= <width>.
    """
    peaks_dsrt = set()
    for idx,(a,b) in enumerate(zip(peaks[:-1], peaks[1:])):
        if b-a < width:
            continue
        s = a
        # print(len(pos_to_rid),b)
        for i in range(a+1,b):
            if len(pos_to_rid[i]) > 0:
                s = i
            if i - s >= width:
                peaks_dsrt.add(idx)
                peaks_dsrt.add(idx+1)
                break
    return peaks_dsrt

def plot(Y_roll, peak_positions, peaks_dyna, peaks_dyna_i, peaks_dyna_o, peaks_dsrt, peaks_vrnc, annotated_brks, pos_to_rid, outpath):
    intervals = [[0,0]]
    margin = 50
    bridge = 200
    for i in annotated_brks:
        Y_roll[i]+= 1
    for idx,i in enumerate(Y_roll):
        if i == 0:
            continue
        if intervals[-1][1]+bridge>=idx:
            intervals[-1][1] = idx
        else:
            intervals[-1][1]+=margin
            intervals.append([idx-margin,idx])
    for i in annotated_brks:
        Y_roll[i]-= 1
    sizes = [max((e-s)/100, 2) for s,e in intervals]
    fig, axes = pp.subplots(1, len(intervals), sharey=True, figsize=(sum(sizes),5), gridspec_kw={'width_ratios': sizes})
    pp.subplots_adjust(left=2/sum(sizes), right=1-2/sum(sizes))
    axes_cov = [ax.twinx() for ax in axes]
    axes_cov[0].get_shared_y_axes().join(*axes)
    fig.text(0.5, 0.04, 'Gene position', ha='center')
    fig.text(0.5, 0.96, '|p|={} |peaks_dyna|={} |peaks_dsrt|={} |peaks_vrnc|={}'.format(len(peak_positions), len(peaks_dyna), len(peaks_dsrt), len(peaks_vrnc)), ha='center')
    fig.text(1/sum(sizes), 0.5, '# of read splicing events', va='center', rotation='vertical')
    plot_data = list()
    plot_data.append(dict(
        X=[peak_positions[i] for i in set(range(len(peak_positions)))-set(peaks_dyna)-set(peaks_dsrt)-set(peaks_vrnc)],
        label='peaks_othr',
        mark='x',
        color='#66c2a5',
    ))
    plot_data.append(dict(
        X=[peak_positions[i] for i in peaks_dyna],
        label='peaks_dyna',
        mark='o',
        color='#e78ac3',
    ))
    plot_data.append(dict(
        X=[peak_positions[i] for i in peaks_dsrt],
        label='peaks_dsrt',
        mark='>',
        color='brown',
    ))
    plot_data.append(dict(
        X=[peak_positions[i] for i in peaks_vrnc],
        label='peaks_vrnc',
        mark='<',
        color='#66c2a5',
    ))
    for d in plot_data:
        d['Y'] = [Y_roll[p] for p in d['X']]
    for idx,ax in enumerate(axes):
        ax.set_xlim(intervals[idx][0],intervals[idx][1])
        for d in plot_data:
            ax.plot(d['X'], d['Y'], d['mark'], label=d['label'], color=d['color'], zorder=10)
        ax.plot(Y_roll, label='Gaussian filtered', color='#fc8d62', zorder=3)
        ax.vlines(annotated_brks,ymin=0, ymax=max(Y_roll)*1.05, colors='#8da0cb', linestyles='dashed', alpha=0.4, linewidth=1, label='Annotation splice points',zorder=4)
        for idxx,(i,j,k) in enumerate(peaks_dyna_o):
            if not intervals[idx][0]<= peak_positions[i] <= intervals[idx][1]:
                continue
            step = idxx%4
            max_h = max(Y_roll)
            height = max_h*0.30+max_h*0.1*step
            ax.plot((peak_positions[i], peak_positions[j]),(height,height),linewidth=3, color='black', alpha=0.6, zorder=0)
            height = max_h*0.40+max_h*0.1*step
            ax.annotate(xy=(peak_positions[i]+(peak_positions[j]-peak_positions[i])*.33,height), s='O: {}'.format(peaks_dyna_o[(i,j,k)]), fontsize=20)
            height = max_h*0.20+max_h*0.1*step
            ax.annotate(xy=(peak_positions[i]+(peak_positions[j]-peak_positions[i])*.33,height), s='I: {}'.format(peaks_dyna_i[(i,j)]), fontsize=20)
        ax.yaxis.tick_left()
        axes_cov[idx].plot([len(p) for p in pos_to_rid], label='Read coverage', color='black', alpha=.4,zorder=2)
        axes_cov[idx].yaxis.tick_right()
    handles, labels   = axes[-1].get_legend_handles_labels()
    handles2, labels2 = axes_cov[-1].get_legend_handles_labels()
    axes[-1].legend(handles+handles2, labels+labels2,bbox_to_anchor=(1+1.5*2/sizes[-1], 0.75), loc='upper left', borderaxespad=0.)
    print(peaks_dyna_i)
    print(peaks_dyna_o)
    fig.savefig(fname=outpath)

def optimize(peaks, C, start, end, low, high, min_read_support_outside):
    cov_mem = dict()
    yea_mem = dict()
    nay_mem = dict()
    amb_mem = dict()

    print('Precomputing coverage mems for {}...'.format((start,end)))
    for i in range(start, end):
        for j in range(i, end+1):
            cov_mem[(i,j)] = (C[j]-C[i])/(peaks[j]-peaks[i]+1)
            yea_mem[(i,j)] = cov_mem[(i,j)] > high
            nay_mem[(i,j)] = cov_mem[(i,j)] < low
            amb_mem[(i,j)] = np.logical_not(np.logical_or(yea_mem[(i,j)],nay_mem[(i,j)]))

    in_mem = dict()
    def inside(i, j):
        if not (i,j) in in_mem:
            if i==j:
                in_mem[(i,j)] = 0
            else:
                in_mem[(i,j)] = -1*amb_mem[(i,j)].sum()
        return in_mem[(i,j)]
    out_mem = dict()
    def outside(i, j, k):
        if not (i,j,k) in out_mem:
            if i==j or j==k:
                out_mem[(i,j,k)] = 0
            else:
                out_mem[(i,j,k)] = sum(np.logical_or(
                    np.logical_and(
                        yea_mem[(i,j)],
                        nay_mem[(j,k)]
                    ),
                    np.logical_and(
                        nay_mem[(i,j)],
                        yea_mem[(j,k)]
                    ),
                ))
                if out_mem[(i,j,k)] < min_read_support_outside:
                    out_mem[(i,j,k)] = -float('inf')
        return out_mem[(i,j,k)]
    D = dict()
    B = dict()
    def dp(i,j,k):
        # memoization
        if (i,j,k) in D or (i,j,k) in B:
            assert (i,j,k) in D and (i,j,k) in B
            return D[(i,j,k)]
        # Base case: i<j<k=END: k is at the end so no more segmentation
        if k == end:
            D[(i,j,k)] = inside(i,j) + outside(i,j,k) + inside(j,k)
            B[(i,j,k)] = (-1,-1,-1)
            return D[(i,j,k)]

        max_b = (-1,-1,-1)
        max_d = float('-inf')
        # Does further segmentation give us better score?
        for k_ in range(k+1, end+1):
            cur_b = (j,k,k_)
            cur_d = inside(i,j) + outside(i,j,k) + dp(*cur_b)
            if cur_d > max_d:
                max_d = cur_d
                max_b = cur_b
        D[(i,j,k)] = max_d
        B[(i,j,k)] = max_b
        return D[(i,j,k)]

    print('DP...')
    # Lower bound on score is no segmentation
    max_d = inside(start,end)
    max_b = (-1,-1,-1)
    for j in range(start+1, end):
        for k in range(j+1, end+1):
            if dp(start,j,k) > max_d:
                max_b = (start,j,k)
                max_d = dp(*max_b)
    print(max_b,max_d)
    return D,B,max_d,max_b,in_mem,out_mem

def get_coverage(peaks, pos_to_rid):
    rids_cnt = 0
    for rids in pos_to_rid:
        for rid in rids:
            if rid > rids_cnt:
                rids_cnt = rid
    rids_cnt += 1
    C = np.zeros((len(peaks)+1, rids_cnt), dtype=np.uint32)
    for idx,(cur,nxt) in enumerate(zip(peaks[:-1],peaks[1:]), start=1):
        for pos in range(cur,nxt):
            for rid in pos_to_rid[pos]:
                C[idx,rid]+=1
    for i in range(1,len(C)):
        C[i] += C[i-1]
    return C

def get_unaligned_gaps(data,brks,t_len,rid_to_intervals,rid_to_len):
    rid_to_unaln_gaps = dict()
    for rid,intervals in rid_to_intervals.items():
        rid_to_unaln_gaps[rid]=list()
        int_idx = 0
        seg_idx = 0
        # print('---',rid, rid_to_len[rid])
        # print(''.join((str(x)for x in data[rid])))
        # print('-'.join((str(x)for x in intervals)))
        # print('...')
        while seg_idx < len(data[rid]):
            eo_idx = -1
            so_idx = len(data[rid])
            # Find the next eo_idx such that data[rid][eo_idx]==1 and data[rid][eo_idx+1]!=1
            while seg_idx < len(data[rid]):
                if data[rid][seg_idx] != 1:
                    break
                eo_idx = seg_idx
                seg_idx += 1
            if eo_idx+1==len(data[rid]):
                continue
            assert data[rid][eo_idx+1]!=1 and (eo_idx==-1 or data[rid][eo_idx]==1)
            # Find the next so_idx such that data[rid][so_idx]==1 and data[rid][so_idx-1]!=1
            while seg_idx < len(data[rid]):
                if data[rid][seg_idx] == 1:
                    so_idx = seg_idx
                    break
                seg_idx += 1
            assert data[rid][so_idx-1]!=1 and (so_idx==len(data[rid]) or data[rid][so_idx]==1)
            gap_ts = brks[eo_idx+1]
            gap_te = brks[so_idx]
            # if gap_ts == 0 or gap_te==t_len:
            #     continue
            # print('(eo_idx, so_idx)',(eo_idx, so_idx))
            # print('(gap_ts, gap_te)',(gap_ts, gap_te))
            gap_qs = 0
            gap_qe = rid_to_len[rid]
            for i in range(int_idx,len(intervals)):
                ((ts,te),(qs,qe)) = intervals[i]
                # print('gqs:',i,intervals[i])
                if ts < gap_ts:
                    gap_qs = qe + (gap_ts-te)
                    int_idx=i
                if ts >= gap_ts:
                    break
            for i in range(int_idx,len(intervals)):
                ((ts,te),(qs,qe)) = intervals[i]
                int_idx=i
                # print('gqe:',i,intervals[i])
                if te >= gap_te:
                    gap_qe = qs - (ts-gap_te)
                    break
            rid_to_unaln_gaps[rid].append((eo_idx+1,so_idx-1,max(0,gap_qe-gap_qs)))
    return rid_to_unaln_gaps

def main():
    args = parse_args()
    assert(args.max_peak_cnt_per_segment>0)
    pos_to_rid,rid_to_len,rid_to_intervals,read_name_to_id,t_len = read_paf(args.paf)
    annotated_brks = annotated_breakpoints(args.tsv)
    if len(rid_to_intervals)>0:
        gene_len = int(open(args.paf).readline().split('\t')[6])
    else:
        gene_len = 1
    X, Y_i, Y_o, Y_a = get_splice(gene_len=gene_len, rid_to_intervals=rid_to_intervals)

    if not args.sigma > 0.0:
        print('Sigma is not positive float. Not applying any gaussian filter')
        Y_roll=Y_a
    else:
        Y_roll=gaussian_filter1d(Y_a,args.sigma)
    peak_positions, _ = find_peaks(Y_roll)
    peak_positions = sorted(set(peak_positions) | {0,t_len})
    print(peak_positions)

    peaks_dsrt = get_desert_bound_peaks(peaks=peak_positions, pos_to_rid=pos_to_rid)
    peaks_vrnc = get_high_var_peaks(peaks=peak_positions, Y_roll=Y_roll, variance_factor=args.variance_factor)
    peaks_fixd = sorted(peaks_dsrt|peaks_vrnc)
    for s,e in zip(peaks_fixd[:-1],peaks_fixd[1:]):
        peak_cnt = (e-s)
        extra_peak_cnt = peak_cnt//args.max_peak_cnt_per_segment
        for i in range(extra_peak_cnt):
            peaks_fixd.append(s+(i+1)*peak_cnt//extra_peak_cnt)
    peaks_fixd = sorted(set(peaks_fixd))
    peaks_dyna = set()
    peaks_dyna_i = dict()
    peaks_dyna_o = dict()

    coverage = get_coverage(peaks=peak_positions, pos_to_rid=pos_to_rid)

    optimizing_args = list()
    for start,end in zip(peaks_fixd[:-1],peaks_fixd[1:]):
        if end-start+1 < 3:
            continue
        print((start,end))
        optimizing_args.append((
            peak_positions,
            coverage,
            start,
            end,
            args.low_threshold,
            args.high_threshold,
            args.min_read_support_outside,
        ))
    with Pool(args.threads) as p:
        for D,B,max_d,max_b,in_mem,out_mem in p.starmap(optimize, optimizing_args):
            while max_b != (-1,-1,-1):
                peaks_dyna.update(max_b)
                peaks_dyna_i[(max_b[0],max_b[1])] = in_mem[(max_b[0],max_b[1])]
                peaks_dyna_o[max_b] = out_mem[max_b]
                print('B:', max_b)
                print('I:', peaks_dyna_i[(max_b[0],max_b[1])])
                print('O:', peaks_dyna_o[max_b])
                max_b = B[max_b]

    peaks_dyna -= set(peaks_fixd)
    peaks_dyna = sorted(peaks_dyna)
    peaks_fina = sorted(set(peaks_dyna) | set(peaks_dsrt) | set(peaks_vrnc) | {0,len(peak_positions)-1})
    plot(
        Y_roll=Y_roll,
        peak_positions=peak_positions,
        peaks_dyna=peaks_dyna,
        peaks_dyna_i=peaks_dyna_i,
        peaks_dyna_o=peaks_dyna_o,
        peaks_dsrt=peaks_dsrt,
        peaks_vrnc=peaks_vrnc,
        annotated_brks=annotated_brks,
        pos_to_rid=pos_to_rid,
        outpath='{}.pdf'.format(args.out_prefix)
    )
    out_file = open('{}.txt'.format(args.out_prefix), 'w+')
    for i in peaks_fina:
        print(peak_positions[i], file=out_file)
    out_file.close()

    out_file = open('{}.data'.format(args.out_prefix), 'w+')
    data = np.array([(coverage[j]-coverage[i])/(peak_positions[j]-peak_positions[i]+1) for i,j in zip(peaks_fina[:-1],peaks_fina[1:])])
    data[data > args.high_threshold] = 1
    data[data < args.low_threshold] = 0
    data[(args.high_threshold > data) & (data > args.low_threshold)] = 2

    data = data.transpose().astype(np.int8)
    for l in data:
        print(''.join([str(int(x)) for x in l]),file=out_file)
    out_file.close()
    out_file = open('{}.gaps'.format(args.out_prefix), 'w+')
    rid_to_unaln_gaps = get_unaligned_gaps(data=data,brks=peak_positions,t_len=t_len,rid_to_intervals=rid_to_intervals,rid_to_len=rid_to_len)
    for rid in range(len(data)):
        print('\t'.join(['{}-{}-{}'.format(*x) for x in rid_to_unaln_gaps[rid]]),file=out_file)
    out_file.close()

if __name__ == "__main__":
    main()
