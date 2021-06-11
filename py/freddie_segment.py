#!/usr/bin/env python3
from multiprocessing import Pool
import os
import glob

from itertools import chain, groupby
import argparse
import re
from math import exp, ceil

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

cinterval_re = '([0-9]+)-([0-9]+)'
cigar_op_re = '([0-9]+)([MIDNSHPX=])'
cigar_re = '((%(c)s)+)' % {'c': cigar_op_re}
rinterval_re = '%(i)s:%(i)s:%(c)s' % {
    'i': '([0-9]+)-([0-9]+)', 'c': cigar_re}
tint_prog = re.compile(r'#%(chr_re)s\t%(cid_re)s\t%(intervals_re)s\t%(read_count_re)s\n$' % {
    'chr_re': '(?P<chr>[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)',
    'cid_re': '(?P<cid>[0-9]+)',
    'intervals_re': '(?P<intervals>%(i)s(,%(i)s)*)' % {'i': cinterval_re},
    'read_count_re': '(?P<read_count>[0-9]+)',
})
read_prog = re.compile(r'%(rid_re)s\t%(name_re)s\t%(chr_re)s\t%(strand_re)s\t%(cid_re)s\t%(intervals_re)s\n$' % {
    'rid_re': '(?P<rid>[0-9]+)',
    'name_re': '(?P<name>[!-?A-~]{1,254})',
    'chr_re': '(?P<chr>[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)',
    'strand_re': '(?P<strand>[+-])',
    'cid_re': '(?P<cid>[0-9]+)',
    'intervals_re': r'(?P<intervals>%(i)s(\t%(i)s)*)' % {'i': rinterval_re},
})
cinterval_prog = re.compile(cinterval_re)
rinterval_prog = re.compile(rinterval_re)
cigar_prog = re.compile(cigar_op_re)

rev_comp = dict(
    A='T',
    C='G',
    G='C',
    T='A'
)
dna_id = dict(
    A='A',
    C='C',
    G='G',
    T='T'
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    parser.add_argument("-s",
                        "--split-dir",
                        type=str,
                        required=True,
                        help="Path to Freddie split directory of the reads")
    parser.add_argument("-o",
                        "--outdir",
                        type=str,
                        default='freddie_segment/',
                        help="Path to output directory. Default: freddie_segment/")
    parser.add_argument("-t",
                        "--threads",
                        type=int,
                        default=1,
                        help="Number of threads for multiprocessing. Default: 1")
    parser.add_argument("-sd",
                        "--sigma",
                        type=float,
                        default=5.0,
                        help="Sigma value for gaussian_filter1d")
    parser.add_argument("-tp",
                        "--threshold-rate",
                        type=float,
                        default=0.90,
                        help="Threshold rate above which the read will be considered as covering a segment. Low threshold is 1-threshold_rate. Anything in between is considered ambigious. Default: 0.9. Note: the stricter threshold for a given segment length will be used.")
    parser.add_argument("-vf",
                        "--variance-factor",
                        type=float,
                        default=3.0,
                        help="The stdev factor to fix a candidate peak. The threshold is set as > mean(non-zero support for splicing postions)+variance_factor*stdev(non-zero support for splicing postions). Default 3.0")
    parser.add_argument("-mps",
                        "--max-problem-size",
                        type=int,
                        default=50,
                        help="Maximum number of candidate breakpoints allowed per segmentation problem")
    parser.add_argument("-lo",
                        "--min-read-support-outside",
                        type=int,
                        default=3,
                        help="Minimum reads support for splice site to support a breakpoint")
    args = parser.parse_args()

    assert(1 >= args.threshold_rate >= 0.5)
    assert(10 > args.variance_factor > 0)
    assert(50 >= args.sigma > 0)
    assert(args.max_problem_size > 3)
    assert(args.min_read_support_outside >= 0)
    assert(args.threads > 0)
    return args


def read_split(split_tsv):
    tints = dict()

    for line in open(split_tsv):
        if line[0] == '#':
            re_dict = tint_prog.match(line).groupdict()
            tint = dict(
                id=int(re_dict['cid']),
                chr=re_dict['chr'],
                intervals=[(int(x[0]), int(x[1]))
                           for x in cinterval_prog.findall(re_dict['intervals'])],
                read_count=int(re_dict['read_count']),
                reads=list()
            )
            assert not tint['id'] in tints, 'Transcriptional interval with id {} is repeated!'.format(
                tint['id'])
            assert all(a[1] < b[0] for a, b in zip(
                tint['intervals'][:-1], tint['intervals'][1:])), (tint['intervals'])
            assert all(s < e for s, e in tint['intervals'])
            tints[tint['id']] = tint
        else:
            re_dict = read_prog.match(line).groupdict()
            read = dict(
                id=int(re_dict['rid']),
                name=re_dict['name'],
                chr=re_dict['chr'],
                strand=re_dict['strand'],
                tint=int(re_dict['cid']),
                intervals=[(
                    int(x[0]),
                    int(x[1]),
                    int(x[2]),
                    int(x[3]),
                    [(int(c[0]), c[1]) for c in cigar_prog.findall(x[4])],
                ) for x in rinterval_prog.findall(re_dict['intervals'])],
            )
            assert all(aet <= bst and aer <= bsr for (_, aet, _, aer, _), (bst, _, bsr, _, _) in zip(
                read['intervals'][:-1], read['intervals'][1:]))
            assert all(st < et and sr < er for (
                st, et, sr, er, _) in read['intervals'])
            tints[read['tint']]['reads'].append(read)
    for tint in tints.values():
        assert len(tint['reads']) == tint['read_count']
    return list(tints.values())


def read_sequence(tint, reads_tsv):
    rid_to_seq = dict()
    for line in open(reads_tsv):
        line = line.rstrip().split('\t')
        rid = int(line[0])
        seq = line[3]
        rid_to_seq[rid] = seq
    assert len(rid_to_seq) == len(tint['reads']), tint
    for read in tint['reads']:
        read['seq'] = rid_to_seq[read['id']]
        read['length'] = len(read['seq'])
    return


def get_cumulative_coverage(candidate_y_idxs, y_idx_to_r_idxs):
    C = np.zeros((len(candidate_y_idxs)+1,
                  y_idx_to_r_idxs.shape[1]), dtype=np.uint32)
    for C_idx, (cur_y_idx, nxt_y_idx) in enumerate(zip(candidate_y_idxs[:-1], candidate_y_idxs[1:]), start=1):
        C[C_idx] = y_idx_to_r_idxs[cur_y_idx:nxt_y_idx].sum(axis=0)
    for C_idx in range(1, len(C)):
        C[C_idx] += C[C_idx-1]
    return C


def refine_segmentation(y_raw, y_idxs, sigma, skip=20, min_internal_splice=20):
    refine_y_idxs = list()
    for s_yidx, e_yidx in zip(y_idxs[:-1], y_idxs[1:]):
        if e_yidx-s_yidx <= 2*skip:
            continue
        i_vals = [x for x in y_raw[s_yidx:e_yidx]]
        for i in range(0, skip):
            i_vals[i] = 0.0
            i_vals[-i-1] = 0.0
        if sum(i_vals) < min_internal_splice:
            continue
        i_gauss = gaussian_filter1d(
            i_vals, sigma, mode='constant', cval=0.0, truncate=1.0)
        for i in find_peaks(i_gauss, distance=skip)[0]:
            if sum(i_gauss[int(round(i-sigma)):int(round(i+sigma+1))]) < min_internal_splice:
                continue
            refine_y_idxs.append(i+s_yidx)
    return refine_y_idxs


def get_high_threshold(seg_len, smoothed_threshold, threshold_rate):
    if seg_len < len(smoothed_threshold):
        h = smoothed_threshold[seg_len]
    else:
        h = threshold_rate
    return h


def smooth_threshold(threshold):
    smooth = list()
    while True:
        x = len(smooth)
        y = threshold/(1 + ((threshold-.5)/.5)*exp(-0.05*x))
        if x > 5 and x*(threshold-y) < 0.5:
            break
        smooth.append(round(y, 2))
        assert len(smooth) < 1000
    return smooth


def forward_thread_cigar(cigar, t_goal, t_pos, q_pos):
    assert t_pos <= t_goal
    cig_idx = 0
    while t_pos < t_goal:
        c, t = cigar[cig_idx]
        c = min(c, t_goal-t_pos)
        if t in ['M', 'X', '=']:
            t_pos += c
            q_pos += c
        if t in ['D']:
            t_pos += c
        if t in ['I']:
            q_pos += c
        cig_idx += 1
    assert t_pos == t_goal
    return q_pos


def get_interval_start(start, read):
    '''
    Finds the first position in the read alignment intervals that aligns to start or after.
    Any (negative) offset is reported as slack.
    '''
    for (t_start, t_end, q_start, q_end, cigar) in read['intervals']:
        # We haven't passed the start location yet
        if t_end < start:
            continue
        if start < t_start:
            q_pos = q_start
            slack = start-t_start
        else:
            q_pos = forward_thread_cigar(
                cigar=cigar, t_goal=start, t_pos=t_start, q_pos=q_start)
            slack = 0
        assert slack <= 0, (slack, t_start, start)
        assert q_start <= q_pos <= q_end, (q_start, q_pos, q_end)
        return q_pos, slack
    assert False


def get_interval_end(end, read):
    '''
    Finds the last location on the read that aligns on or before end.
    Any (negative) offset is reported as slack.
    '''
    for (t_start, t_end, q_start, q_end, cigar) in reversed(read['intervals']):
        # We haven't passed end yet
        if t_start > end:
            continue
        # the right most interval that covers end, ends before end
        if t_end < end:
            q_pos = q_end
            slack = t_end-end
        else:
            q_pos = forward_thread_cigar(
                cigar=cigar, t_goal=end, t_pos=t_start, q_pos=q_start)
            slack = 0
        assert slack <= 0, (slack, t_end, end)
        assert 0 <= q_pos <= q_end, (q_start, q_pos, q_end)
        return q_pos, slack
    assert False


def find_longest_poly(seq, s, e, step, match_score=1, mismatch_score=-2, char='A'):
    if e-s == 0:
        return
    if seq[s] == char:
        scores = [match_score]
    else:
        scores = [0]
    for m in (match_score if c == char else mismatch_score for c in seq[s+step:e:step]):
        scores.append(max(0, scores[-1]+m))
    for k, g in groupby(enumerate(scores), lambda x: x[1] > 0):
        if not k:
            continue
        i, S = list(zip(*g))
        max_s, max_i = max(zip(S, i))
        l = max_i+1-i[0]
        yield i[0], l, seq[s:e:step][i[0]:i[0]+l].count(char)/l


def get_unaligned_gaps_and_polyA(read, segs):
    read['gaps'] = set()
    if not 1 in read['data']:
        return
    intervals = list()
    for d, group in groupby(enumerate(read['data']), lambda x: x[1]):
        if d != 1:
            continue
        group = list(group)
        f_seg_idx = group[0][0]
        l_seg_idx = group[-1][0]
        intervals.append((f_seg_idx, l_seg_idx))
    assert len(intervals) > 0, read['data']
    (f_seg_idx, _) = intervals[0]
    start = segs[f_seg_idx][0]
    q_ssc_pos, _ = get_interval_start(start=start, read=read)
    (_, l_seg_idx) = intervals[-1]
    end = segs[l_seg_idx][1]
    q_esc_pos, _ = get_interval_end(end=end, read=read)
    assert 0 <= q_ssc_pos <= q_esc_pos <= read['length'], (
        q_ssc_pos, q_esc_pos, read['length'], start, end, segs, read)
    s_polys = list()
    for char in ['A', 'T']:
        s = 0
        e = q_ssc_pos
        step = 1
        sc_char = char
        if read['strand'] == '-':
            s = -s -1
            e = -e -1
            step = -1
            sc_char = rev_comp[char]
        for i, l, p in find_longest_poly(read['seq'], s=s, e=e, step=step, char=sc_char):
            if l < 20 or p < 0.85:
                continue
            assert 0 <= i < q_ssc_pos, (i, q_ssc_pos, read['length'])
            s_polys.append((i, l, p, char))
    if len(s_polys) > 0:
        i, l, p, char = max(s_polys, key=lambda x: x[2])
        poly_to_gene_gap_size = q_ssc_pos-i-l
        assert 0 <= poly_to_gene_gap_size < q_ssc_pos
        read['gaps'].add(
            'S{}_{}:{}'.format(char, l, poly_to_gene_gap_size)
        )
        read['gaps'].add(
            'SSC:{}'.format(i)
        )
    else:
        read['gaps'].add(
            'SSC:{}'.format(q_ssc_pos)
        )
    e_polys = list()
    for char in ['A', 'T']:
        s = q_esc_pos
        e = read['length']
        step = 1
        sc_char = char
        if read['strand'] == '-':
            s = -s -1
            e = -e -1
            step = -1
            sc_char = rev_comp[char]
        for i, l, p in find_longest_poly(read['seq'], s=s, e=e, step=step, char=sc_char):
            if l < 20 or p < 0.85:
                continue
            assert 0 <= i < read['length'] - \
                q_esc_pos, (i, q_esc_pos, read['length'])
            e_polys.append((i, l, p, char))
    if len(e_polys) > 0:
        i, l, p, char = max(e_polys, key=lambda x: x[2])
        poly_to_gene_gap_size = i
        assert 0 <= poly_to_gene_gap_size < read['length'] - \
            q_esc_pos, (q_esc_pos, i, l, p,
                        read['length'], poly_to_gene_gap_size)
        read['gaps'].add(
            'E{}_{}:{}'.format(char, l, poly_to_gene_gap_size)
        )
        read['gaps'].add(
            'ESC:{}'.format(read['length']-q_esc_pos-poly_to_gene_gap_size)
        )
        assert read['length']-q_esc_pos-poly_to_gene_gap_size > 0
    else:
        read['gaps'].add(
            'ESC:{}'.format(read['length']-q_esc_pos)
        )
    for i1, i2 in zip(intervals[:-1], intervals[1:]):
        (_, i1_l_seg_idx) = i1
        i1_end = segs[i1_l_seg_idx][1]
        q_gap_start, start_slack = get_interval_end(end=i1_end, read=read)
        (i2_f_seg_idx, _) = i2
        i2_start = segs[i2_f_seg_idx][0]
        q_gap_end, end_slack = get_interval_start(start=i2_start, read=read)
        assert 0 < q_gap_start <= q_gap_end < read['length'], (
            q_gap_start, q_gap_end, read['length'])
        q_gap_size = q_gap_end-q_gap_start
        q_gap_size = max(0, q_gap_size+start_slack+end_slack)
        assert 0 <= q_gap_size < read['length'], (
            q_gap_size, start_slack, end_slack)
        assert i1_l_seg_idx < i2_f_seg_idx
        read['gaps'].add(
            '{}-{}:{}'.format(i1_l_seg_idx, i2_f_seg_idx, q_gap_size),
        )
    read['gaps'] = sorted(read['gaps'])


def optimize(candidate_y_idxs, C, start, end, smoothed_threshold, threshold_rate, read_support):
    cov_mem = dict()
    yea_mem = dict()
    nay_mem = dict()
    amb_mem = dict()

    # print('Precomputing coverage mems for {}...'.format((start,end)))
    for i in range(start, end):
        for j in range(i, end+1):
            seg_len = (candidate_y_idxs[j]-candidate_y_idxs[i]+1)
            cov_mem[(i, j)] = (C[j]-C[i])/seg_len
            h = get_high_threshold(seg_len, smoothed_threshold, threshold_rate)
            l = 1-h
            yea_mem[(i, j)] = cov_mem[(i, j)] > h
            nay_mem[(i, j)] = cov_mem[(i, j)] < l
            amb_mem[(i, j)] = np.logical_not(
                np.logical_or(yea_mem[(i, j)], nay_mem[(i, j)]))

    in_mem = dict()

    def inside(i, j):
        if not (i, j) in in_mem:
            if i == j:
                in_mem[(i, j)] = 0
            else:
                in_mem[(i, j)] = -1*amb_mem[(i, j)].sum()
        return in_mem[(i, j)]
    out_mem = dict()

    def outside(i, j, k):
        if not (i, j, k) in out_mem:
            if i == j or j == k:
                out_mem[(i, j, k)] = 0
            else:
                out_mem[(i, j, k)] = np.sum(np.logical_or(
                    np.logical_and(
                        yea_mem[(i, j)],
                        nay_mem[(j, k)]
                    ),
                    np.logical_and(
                        nay_mem[(i, j)],
                        yea_mem[(j, k)]
                    ),
                ))
                if out_mem[(i, j, k)] < read_support:
                    out_mem[(i, j, k)] = float('-inf')
        return out_mem[(i, j, k)]
    D = dict()
    B = dict()

    def dp(i, j, k):
        # memoization
        if (i, j, k) in D or (i, j, k) in B:
            assert (i, j, k) in D and (i, j, k) in B
            return D[(i, j, k)]
        max_b = (-1, -1, -1)
        max_d = float('-inf')
        # Segment too small: y_idx[j]-y_idx[i] < 5 or
        if candidate_y_idxs[j]-candidate_y_idxs[i] < 5 or candidate_y_idxs[k]-candidate_y_idxs[j] < 5:
            D[(i, j, k)] = max_d
            B[(i, j, k)] = max_b
            return D[(i, j, k)]
        # Base case: i<j<k=END: k is at the end so no more segmentation
        if k == end:
            D[(i, j, k)] = inside(i, j) + outside(i, j, k) + inside(j, k)
            B[(i, j, k)] = (-1, -1, -1)
            return D[(i, j, k)]
        # Does further segmentation give us better score?
        for k_ in range(k+1, end+1):
            cur_b = (j, k, k_)
            cur_d = inside(i, j) + outside(i, j, k) + dp(*cur_b)
            if cur_d > max_d:
                max_d = cur_d
                max_b = cur_b
        D[(i, j, k)] = max_d
        B[(i, j, k)] = max_b
        return D[(i, j, k)]

    # print('DP...')
    # Lower bound on score is no segmentation
    max_d = inside(start, end)
    max_b = (-1, -1, -1)
    for j in range(start+1, end):
        for k in range(j+1, end+1):
            if dp(start, j, k) > max_d:
                max_b = (start, j, k)
                max_d = dp(*max_b)
    # print(max_b,max_d)
    return D, B, max_d, max_b, in_mem, out_mem


def run_optimize(candidate_y_idxs, fixed_c_idxs, coverage, smoothed_threshold, threshold_rate, min_read_support_outside):
    final_c_idxs = set(fixed_c_idxs)
    for idx, (start, end) in enumerate(zip(fixed_c_idxs[:-1], fixed_c_idxs[1:])):
        # print('{}/{}: {}'.format(idx+1, len(fixed_c_idxs)-1, end-start+1))
        D, B, max_d, max_b, in_mem, out_mem = optimize(
            candidate_y_idxs=candidate_y_idxs,
            C=coverage,
            start=start,
            end=end,
            smoothed_threshold=smoothed_threshold,
            threshold_rate=threshold_rate,
            # threshold_abs            = threshold_abs,
            read_support=min_read_support_outside,
        )
        while max_b != (-1, -1, -1):
            final_c_idxs.update(max_b)
            # print('D:', D[max_b])
            # print('B:', max_b)
            # print('I:', in_mem[(max_b[0],max_b[1])])
            # print('O:', out_mem[max_b])
            max_b = B[max_b]
    return sorted(final_c_idxs)


def non_desert(y, jump=10):
    l = list()
    for k, group in groupby(enumerate(y), lambda x: x[1] > 0):
        if not k:
            continue
        group = list(group)
        f_idx = group[0][0]
        l_idx = group[-1][0]
        if len(l) == 0:
            l.append([f_idx, l_idx])
        elif l_idx - l[-1][-1] < jump:
            l[-1][-1] = l_idx
        else:
            l.append([f_idx, l_idx])
    return l


def candidates_from_peaks(y, f_y_idx, l_y_idx):
    # print('f_y_idx, l_y_idx',f_y_idx, l_y_idx)
    c, _ = find_peaks(y[f_y_idx:l_y_idx+1])
    c = [f_y_idx+y_idx for y_idx in c]
    return c


def candidates_from_window(y, f_y_idx, l_y_idx, window=5):
    # print('f_y_idx, l_y_idx',f_y_idx, l_y_idx)
    c = list()
    for y_idx in range(f_y_idx, l_y_idx+1, window):
        max_y_idx = y_idx + np.argmax(y[y_idx:y_idx+window])
        if y[max_y_idx] > 0.1:
            c.append(max_y_idx)
        peaks = find_peaks(y[y_idx:y_idx+window], distance=window)[0]
        assert len(peaks) <= 1
        if len(peaks) > 0:
            peak_y_idx = peaks[0] + y_idx
            if peak_y_idx != max_y_idx and y[peak_y_idx] > 0.1:
                c.append(peak_y_idx)
    c.sort()
    c_keep = [True for _ in c]
    c_idx = 0
    # print('c', c)
    while c_idx < len(c)-1:
        if c[c_idx]+1 == c[c_idx+1]:
            if y[c[c_idx]] > y[c[c_idx+1]]:
                c_keep[c_idx+1] = False
            else:
                c_keep[c_idx] = False
            c_idx += 1
        c_idx += 1
    c = [c[c_idx] for c_idx in range(len(c)) if c_keep[c_idx]]
    return c


def break_large_problems(candidate_y_idxs, fixed_c_idxs, y, max_problem_size, window=5):
    fixed_c_idxs_pairs = sorted(fixed_c_idxs)
    fixed_c_idxs_pairs = [(s, e) for s, e in zip(
        fixed_c_idxs_pairs[:-1], fixed_c_idxs_pairs[1:])]
    new_problems_total_count = 0
    for c_idx_s, c_idx_e in fixed_c_idxs_pairs:
        problem_size = c_idx_e-c_idx_s+1
        if problem_size <= max_problem_size:
            continue
        new_problems_count = ceil(problem_size / max_problem_size)
        new_problems_size = problem_size/new_problems_count
        new_problems_total_count += new_problems_count-1
        for i in range(1, new_problems_count):
            mid_anchor = int(c_idx_s+i*new_problems_size)
            max_c_idx_y_v = float('-inf')
            max_c_idx = None
            for c_idx in range(mid_anchor-window, mid_anchor+window):
                if y[candidate_y_idxs[c_idx]] > max_c_idx_y_v:
                    max_c_idx_y_v = y[candidate_y_idxs[c_idx]]
                    max_c_idx = c_idx
            assert max_c_idx_y_v > 0
            fixed_c_idxs.add(max_c_idx)
    return fixed_c_idxs, new_problems_total_count


def process_splicing_data(tint):
    pos_to_Yy_idx = dict()
    Yy_idx_to_pos = list()
    Yy_idx_to_r_idxs = list()
    Y_raw = list()
    # print('Generating index for tint {} with {} intervals'.format(tint['id'], tint['intervals']))
    for s, e in tint['intervals']:
        y_idx_to_pos = list()
        for p in range(s, e+1):
            assert not p in pos_to_Yy_idx
            pos_to_Yy_idx[p] = (len(Yy_idx_to_pos), len(y_idx_to_pos))
            y_idx_to_pos.append(p)
        Yy_idx_to_pos.append(y_idx_to_pos)
        Yy_idx_to_r_idxs.append(
            np.zeros((len(y_idx_to_pos), len(tint['reads'])), dtype=bool))
        Y_raw.append(np.zeros(len(y_idx_to_pos)))
    assert len(pos_to_Yy_idx) == sum(len(y_idx_to_pos)
                                     for y_idx_to_pos in Y_raw)
    # print('Building per read coverage data for tint {}'.format(tint['id']))
    for r_idx, read in enumerate(tint['reads']):
        read['data'] = list()
        for ts, te, _, _, cigar in read['intervals']:
            Y_idx_s, y_idx_s = pos_to_Yy_idx[ts]
            Y_idx_e, y_idx_e = pos_to_Yy_idx[te]
            assert Y_idx_s == Y_idx_e, (Y_idx_s,Y_idx_e)
            Y_idx = Y_idx_s
            Y_raw[Y_idx][y_idx_s] += 1
            Y_raw[Y_idx][y_idx_e] += 1
            Yy_idx_to_r_idxs[Y_idx][y_idx_s:y_idx_e,r_idx] = True
    return (
        Yy_idx_to_pos,
        Yy_idx_to_r_idxs,
        Y_raw,
    )


def run_segment(segment_args):
    (
        split_dir,
        outdir,
        contig,
        tint_id,
        sigma,
        smoothed_threshold,
        threshold_rate,
        variance_factor,
        max_problem_size,
        min_read_support_outside
    ) = segment_args
    tints = read_split(
        split_tsv='{}/{}/split_{}_{}.tsv'.format(split_dir, contig, contig, tint_id))
    assert len(tints) == 1
    tint = tints[0]
    read_sequence(
        tint=tint, reads_tsv='{}/{}/reads_{}_{}.tsv'.format(split_dir, contig, contig, tint_id))
    segment(
        tint,
        sigma,
        smoothed_threshold,
        threshold_rate,
        variance_factor,
        max_problem_size,
        min_read_support_outside
    )

    out_file = open('{}/{}/segment_{}_{}.tsv'.format(outdir,
                                                     contig, contig, tint_id), 'w+')
    record = list()
    record.append('#{}'.format(tint['chr']))
    record.append(str(tint['id']))
    record.append(','.join(map(str, tint['final_positions'])))
    out_file.write('\t'.join(record))
    out_file.write('\n')
    for read in tint['reads']:
        record = list()
        record.append(str(read['id']))
        record.append(read['name'])
        record.append(read['chr'])
        record.append(read['strand'])
        record.append(str(read['tint']))
        record.append(''.join(map(str, read['data'])))
        record.append(''.join('{},'.format(g) for g in read['gaps']))
        out_file.write('\t'.join(record))
        out_file.write('\n')
    out_file.close()


def segment(tint, sigma, smoothed_threshold, threshold_rate, variance_factor, max_problem_size, min_read_support_outside):
    (
        Yy_idx_to_pos,
        Yy_idx_to_r_idxs,
        Y_raw,
    ) = process_splicing_data(tint)

    Y = [gaussian_filter1d(y, sigma, truncate=4.0) for y in Y_raw]
    Y_none_zero_vals = np.array([v for y in Y for v in y if v > 0])
    variance_threhsold = Y_none_zero_vals.mean() + variance_factor * \
        Y_none_zero_vals.std()

    tint['final_positions'] = list()
    for Y_idx, y in enumerate(Y):
        candidate_y_idxs = list()
        fixed_c_idxs = set()
        for f_y_idx, l_y_idx in [[0, len(y)-1]]:
            fixed_c_idxs.add(len(candidate_y_idxs))
            current_candidates = candidates_from_peaks(y, f_y_idx, l_y_idx)
            current_candidates = sorted(
                set(current_candidates) | {f_y_idx, l_y_idx})
            candidate_y_idxs.extend(current_candidates)
            fixed_c_idxs.add(len(candidate_y_idxs)-1)
        for c_idx, y_idx in enumerate(candidate_y_idxs):
            if y[y_idx] > variance_threhsold:
                fixed_c_idxs.add(c_idx)
        # fixed_c_idxs, new_problems_total_count = break_large_problems(
        #     candidate_y_idxs, fixed_c_idxs, y, max_problem_size)
        fixed_c_idxs = sorted(fixed_c_idxs)
        for s, e in zip(fixed_c_idxs[:-1], fixed_c_idxs[1:]):
            assert e-s+1 <= max_problem_size+5

        cumulative_coverage = get_cumulative_coverage(
            candidate_y_idxs, Yy_idx_to_r_idxs[Y_idx])
        # print('Optimizing tint {} with:\n\tcandidate {} locations: {}\n\tfixed {} loations: {}'.format(tint['id'], len(candidate_y_idxs),candidate_y_idxs,len(fixed_c_idxs),[candidate_y_idxs[c_idx] for c_idx in fixed_c_idxs]))
        final_c_idxs = run_optimize(
            candidate_y_idxs=candidate_y_idxs,
            fixed_c_idxs=fixed_c_idxs,
            coverage=cumulative_coverage,
            smoothed_threshold=smoothed_threshold,
            threshold_rate=threshold_rate,
            min_read_support_outside=min_read_support_outside,
        )
        final_y_idxs = [candidate_y_idxs[c_idx] for c_idx in final_c_idxs]
        refine_y_idxs = refine_segmentation(Y_raw[Y_idx], final_y_idxs, sigma)
        final_y_idxs.extend(refine_y_idxs)
        final_y_idxs.sort()
        tint['final_positions'].extend(
            [Yy_idx_to_pos[Y_idx][y_idx] for y_idx in final_y_idxs])
        cumulative_coverage = get_cumulative_coverage(
            final_y_idxs, Yy_idx_to_r_idxs[Y_idx])
        for seg_idx, (s_yidx, e_yidx) in enumerate(zip(final_y_idxs[:-1], final_y_idxs[1:])):
            seg_len = Yy_idx_to_pos[Y_idx][e_yidx] - \
                Yy_idx_to_pos[Y_idx][s_yidx]+1
            h = get_high_threshold(seg_len, smoothed_threshold, threshold_rate)
            l = 1-h
            assert seg_len == e_yidx-s_yidx+1
            for r_idx, read in enumerate(tint['reads']):
                cov_ratio = (
                    cumulative_coverage[seg_idx+1][r_idx]-cumulative_coverage[seg_idx][r_idx])/seg_len
                assert 0 <= cov_ratio <= 1, (r_idx, seg_idx, s_yidx, e_yidx, seg_len, (
                    cumulative_coverage[seg_idx][r_idx]-cumulative_coverage[seg_idx+1][r_idx]), cumulative_coverage[seg_idx][r_idx], cumulative_coverage[seg_idx+1][r_idx], cov_ratio, read)
                if cov_ratio > h:
                    read['data'].append(1)
                elif cov_ratio < l:
                    read['data'].append(0)
                else:
                    read['data'].append(2)
        for r_idx, read in enumerate(tint['reads']):
            read['data'].append(0)
    tint['segs'] = [(s, e) for s, e in zip(
        tint['final_positions'][:-1], tint['final_positions'][1:])]
    # print('Extracting unaligned gaps and polyA tail data from reads for tint {}'.format(tint['id']))
    for read in tint['reads']:
        read['data'].pop()
        assert len(read['data']) == len(
            tint['segs']), (read['data'], tint['segs'])
        get_unaligned_gaps_and_polyA(read=read, segs=tint['segs'])
    return tint['id']


def main():
    args = parse_args()
    args.split_dir = args.split_dir.rstrip('/')

    segment_args = list()
    for contig in os.listdir(args.split_dir):
        if not os.path.isdir('{}/{}'.format(args.split_dir, contig)):
            continue
        os.makedirs('{}/{}'.format(args.outdir, contig), exist_ok=False)
        for tint_id in glob.iglob('{}/{}/split_*.tsv'.format(args.split_dir, contig)):
            tint_id = int(tint_id[:-4].split('/')[-1].split('_')[-1])
            segment_args.append((
                args.split_dir,
                args.outdir,
                contig,
                tint_id,
                args.sigma,
                smooth_threshold(threshold=args.threshold_rate),
                args.threshold_rate,
                args.variance_factor,
                args.max_problem_size,
                args.min_read_support_outside,
            ))
    if args.threads > 1:
        p = Pool(args.threads)
        for idx, tint_id in enumerate(p.imap_unordered(run_segment, segment_args, chunksize=10)):
            if not idx % ceil(len(segment_args)/100) == 0:
                continue
            print('[freddie_segment] Done with {}/{} tints ({:.1%})'.format(
                idx,
                len(segment_args),
                idx/len(segment_args)),
            )
    else:
        for idx, tint_id in enumerate(map(run_segment, segment_args)):
            if not idx % ceil(len(segment_args)/100) == 0:
                continue
            print('[freddie_segment] Done with {}/{} tints ({:.1%})'.format(
                idx,
                len(segment_args),
                idx/len(segment_args)),
            )


if __name__ == "__main__":
    # import cProfile
    # import pstats
    # profiler = cProfile.Profile()
    # profiler.enable()
    main()
    # profiler.disable()
    # stats = pstats.Stats(profiler, stream=open(
    #     'segment.no_RC.cprof', 'w+')).sort_stats('tottime')
    # stats.print_stats()
