#!/usr/bin/env python3
from multiprocessing import Pool

from itertools import starmap
import argparse
import re
from math import ceil

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    parser.add_argument("-s",
                        "--split-tsv",
                        type=str,
                        required=True,
                        help="Path to Freddie split TSV file of the reads")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default='freddie_segment.tsv',
                        help="Path to output file. Default: freddie_segment.tsv")
    parser.add_argument("-c",
                        "--threads",
                        type=int,
                        default=1,
                        help="Number of threads for multiprocessing. Default: 1")
    parser.add_argument("-sd",
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
    parser.add_argument("-v",
                        "--variance-factor",
                        type=float,
                        default=3.0,
                        help="The sigma factor to fix a candidate peak. The threshold is set as > mean+3*variance_factor. Default 3.0")
    parser.add_argument("-hp",
                        "--max-candidates-per-seg",
                        type=int,
                        default=50,
                        help="Maximum number of candidate breakpoints allowed per segmentation problem")
    parser.add_argument("-lo",
                        "--min-read-support-outside",
                        type=int,
                        default=3,
                        help="Minimum reads support for splice site to support a breakpoint")
    args = parser.parse_args()


    assert(1 >= args.high_threshold > args.low_threshold >= 0)
    assert(args.variance_factor>0)
    assert(args.sigma>0)
    assert(args.max_candidates_per_seg>0)
    assert(args.min_read_support_outside>=0)
    assert(args.threads>0)
    return args

def read_split(split_tsv):
    tints = list()

    cinterval_re   = '([0-9]+)-([0-9]+)'
    cigar_op_re = '([0-9]+)([MIDNSHPX=])'
    cigar_re = '((%(c)s)+)' % {'c':cigar_op_re}
    rinterval_re   = '%(i)s:%(i)s:%(c)s' % {'i':'([0-9]+)-([0-9]+)', 'c':cigar_re}
    tint_prog   = re.compile(r'#%(chr_re)s\t%(cid_re)s\t%(intervals_re)s\t%(read_count_re)s\n$' % {
        'chr_re'        : '(?P<chr>[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)',
        'cid_re'        : '(?P<cid>[0-9]+)',
        'intervals_re'  : '(?P<intervals>%(i)s(,%(i)s)*)' % {'i':cinterval_re},
        'read_count_re' : '(?P<read_count>[0-9]+)',
    })
    read_prog      = re.compile(r'%(rid_re)s\t%(name_re)s\t%(chr_re)s\t%(strand_re)s\t%(cid_re)s\t%(intervals_re)s\n$' % {
        'rid_re'       : '(?P<rid>[0-9]+)',
        'name_re'      : '(?P<name>[!-?A-~]{1,254})',
        'chr_re'       : '(?P<chr>[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)',
        'strand_re'    : '(?P<strand>[+-])',
        'cid_re'       : '(?P<cid>[0-9]+)',
        'intervals_re' : r'(?P<intervals>%(i)s(\t%(i)s)*)' % {'i' : rinterval_re},
    })
    cinterval_prog  = re.compile(cinterval_re)
    rinterval_prog  = re.compile(rinterval_re)
    cigar_prog     = re.compile(cigar_op_re)
    for line in open(split_tsv):
        if line[0]=='#':
            re_dict = tint_prog.match(line).groupdict()
            tint = dict(
                id         = int(re_dict['cid']),
                chr        = re_dict['chr'],
                intervals  = [(int(x[0]),int(x[1])) for x in cinterval_prog.findall(re_dict['intervals'])],
                read_count = int(re_dict['read_count']),
                reads      = list()
            )
            if tint['id'] >= len(tints):
                tints.extend([None]*(tint['id']-len(tints)+1))
            assert tints[tint['id']] == None, 'Transcriptional interval with id {} is repeated!'.format(tint['id'])
            assert all(a[1]<=b[0] for a,b in zip(tint['intervals'][:-1],tint['intervals'][1:])),tint['intervals']
            assert all(s<e for s,e in tint['intervals'])
            tints[tint['id']] = tint
        else:
            re_dict = read_prog.match(line).groupdict()
            re_dict['rid']        = int(re_dict['rid'])
            re_dict['cid']        = int(re_dict['cid'])
            read = dict(
                id         = int(re_dict['rid']),
                name       = re_dict['name'],
                chr        = re_dict['chr'],
                strand     = re_dict['strand'],
                tint       = int(re_dict['cid']),
                intervals  = [(
                    int(x[0]),
                    int(x[1]),
                    int(x[2]),
                    int(x[3]),
                    [(int(c[0]),c[1]) for c in cigar_prog.findall(x[4])],
                ) for x in rinterval_prog.findall(re_dict['intervals'])],
            )
            assert all(aet<=bst and aer<=bsr for (_,aet,_,aer,_),(bst,_,bsr,_,_) in zip(read['intervals'][:-1],read['intervals'][1:]))
            assert all(st<et and sr<er for (st,et,sr,er,_) in read['intervals'])
            tints[read['tint']]['reads'].append(read)
    for x in tints:
        assert x!=None
        len(x['reads'])==x['read_count']
    return tints

def optimize(candidate_y_idxs, C, start, end, low, high, read_support):
    cov_mem = dict()
    yea_mem = dict()
    nay_mem = dict()
    amb_mem = dict()

    # print('Precomputing coverage mems for {}...'.format((start,end)))
    for i in range(start, end):
        for j in range(i, end+1):
            cov_mem[(i,j)] = (C[j]-C[i])/(candidate_y_idxs[j]-candidate_y_idxs[i]+1)
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
                if out_mem[(i,j,k)] < read_support:
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

    # print('DP...')
    # Lower bound on score is no segmentation
    max_d = inside(start,end)
    max_b = (-1,-1,-1)
    for j in range(start+1, end):
        for k in range(j+1, end+1):
            if dp(start,j,k) > max_d:
                max_b = (start,j,k)
                max_d = dp(*max_b)
    # print(max_b,max_d)
    return D,B,max_d,max_b,in_mem,out_mem


def run_optimize(candidate_y_idxs, fixed_c_idxs, coverage, low_threshold, high_threshold, min_read_support_outside):
    final_y_idxs = set(fixed_c_idxs)
    for start,end in zip(fixed_c_idxs[:-1],fixed_c_idxs[1:]):
        D,B,max_d,max_b,in_mem,out_mem = optimize(
            candidate_y_idxs         = candidate_y_idxs,
            C                        = coverage,
            start                    = start,
            end                      = end,
            low                      = low_threshold,
            high                     = high_threshold,
            read_support             = min_read_support_outside,
        )
        while max_b != (-1,-1,-1):
            final_y_idxs.update(max_b)
            # print('B:', max_b)
            # print('I:', in_mem[(max_b[0],max_b[1])])
            # print('O:', out_mem[max_b])
            max_b = B[max_b]
    return sorted(final_y_idxs)


def get_cumulative_coverage(candidate_y_idxs, y_idx_to_r_idxs):
    C = np.zeros((len(candidate_y_idxs)+1, y_idx_to_r_idxs.shape[1]), dtype=np.uint32)
    for C_idx,(cur_y_idx,nxt_y_idx) in enumerate(zip(candidate_y_idxs[:-1],candidate_y_idxs[1:]), start=1):
        C[C_idx] = y_idx_to_r_idxs[cur_y_idx:nxt_y_idx].sum(axis=0)
    for C_idx in range(1,len(C)):
        C[C_idx] += C[C_idx-1]
    return C

def segment(segment_args):
    tint, sigma, low_threshold, high_threshold, variance_factor, max_candidates_per_seg, min_read_support_outside, threads = segment_args
    pos_to_Yy_idx = dict()
    Yy_idx_to_pos = list()
    Yy_idx_to_r_idxs = list()
    Y = list()
    for s,e in tint['intervals']:
        y_idx_to_pos = list()
        for p in range(s,e):
            pos_to_Yy_idx[p]=(len(Yy_idx_to_pos),len(y_idx_to_pos))
            y_idx_to_pos.append(p)
        Yy_idx_to_pos.append(y_idx_to_pos)
        Yy_idx_to_r_idxs.append(np.zeros((len(y_idx_to_pos), len(tint['reads'])), dtype=bool))
        Y.append(np.zeros(len(y_idx_to_pos)))
    assert len(pos_to_Yy_idx)==sum(len(y_idx_to_pos) for y_idx_to_pos in Y)
    for r_idx,read in enumerate(tint['reads']):
        for ts,te,_,_,_ in read['intervals']:
            Y_idx,y_idx = pos_to_Yy_idx[ts]
            Y[Y_idx][y_idx] += 1
            Y_idx,y_idx = pos_to_Yy_idx[te-1]
            Y[Y_idx][y_idx] += 1
            for pos in range(ts,te):
                Y_idx,y_idx = pos_to_Yy_idx[pos]
                Yy_idx_to_r_idxs[Y_idx][y_idx][r_idx] = True
    Y = [gaussian_filter1d(y,sigma) for y in Y]
    Y_none_zero_vals = np.array([v for y in Y for v in y if v > 0])
    variance_threhsold = Y_none_zero_vals.mean() + variance_factor*Y_none_zero_vals.std()

    run_optimize_args = list()
    for Y_idx,y in enumerate(Y):
        candidate_y_idxs,_ = find_peaks(y)
        candidate_y_idxs = sorted(set(candidate_y_idxs) | {0,len(y)-1})
        fixed_c_idxs = {0,len(candidate_y_idxs)-1}
        for c_idx,y_idx in enumerate(candidate_y_idxs):
            if y[y_idx] > variance_threhsold:
                fixed_c_idxs.add(c_idx)
        for s,e in zip(sorted(fixed_c_idxs)[:-1],sorted(fixed_c_idxs)[1:]):
            candidate_count = (e-s)-1
            if candidate_count <= max_candidates_per_seg:
                continue
            extra = candidate_count//max_candidates_per_seg
            brks = np.round(np.linspace(s, e, extra+2)).astype(int)[1:-1]
            fixed_c_idxs.update(brks)
            # print('Add {} extra breakspoints ({}) between {}'.format(extra, brks,(s,e)))
        fixed_c_idxs = sorted(fixed_c_idxs)

        cumulative_coverage = get_cumulative_coverage(candidate_y_idxs, Yy_idx_to_r_idxs[Y_idx])
        run_optimize_args.append((
            candidate_y_idxs,
            fixed_c_idxs,
            cumulative_coverage,
            low_threshold,
            high_threshold,
            min_read_support_outside
        ))
    final_positions = list()
    for Y_idx,final_y_idxs in enumerate(starmap(run_optimize, run_optimize_args)):
        final_positions.extend([Yy_idx_to_pos[Y_idx][y_idx] for y_idx in final_y_idxs])
    return tint['id'],final_positions

def main():
    args = parse_args()

    tints = read_split(args.split_tsv)
    sub_process_threads = min(1, args.threads)
    sup_process_threads = max(1, args.threads//1)
    assert sub_process_threads*sup_process_threads <= args.threads
    # print(sub_process_threads)
    # print(sup_process_threads)
    segment_args = list()
    for tint in tints:
        segment_args.append((
            tint,
            args.sigma,
            args.low_threshold,
            args.high_threshold,
            args.variance_factor,
            args.max_candidates_per_seg,
            args.min_read_support_outside,
            sub_process_threads,
        ))
    with Pool(sup_process_threads) as p:
        for idx,(tint_idx,final_positions) in enumerate(p.imap_unordered(segment, segment_args, chunksize=20)):
            print('Done with {}-th transcriptional multi-intervals ({}/{})'.format(tint_idx, idx+1,len(tints)))
            tints[idx]['final_positions']=final_positions

    for tint in tints:
        print(tint['final_positions'])

if __name__ == "__main__":
    main()
