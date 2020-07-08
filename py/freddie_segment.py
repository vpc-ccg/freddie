#!/usr/bin/env python3
from multiprocessing import Pool

from itertools import chain,groupby
import argparse
import re
from math import ceil,exp

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
    parser.add_argument("-r",
                        "--reads",
                        nargs="+",
                        type=str,
                        required=True,
                        help="Space separated paths to reads in FASTQ or FASTA format used to extract polyA tail information")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default='freddie_segment.tsv',
                        help="Path to output file. Default: freddie_segment.tsv")
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
    parser.add_argument("-td",
                        "--threshold",
                        type=float,
                        default=0.90,
                        help="Threshold above which the read will be considered as covering a segment. Low threshold is 1-threshold. Anything in between is considered ambigious. Default: 0.9")
    parser.add_argument("-vf",
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

    assert(1 >= args.threshold >= 0.5)
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
            assert all(a[1]<b[0] for a,b in zip(tint['intervals'][:-1],tint['intervals'][1:])),(tint['intervals'])
            assert all(s<e for s,e in tint['intervals'])
            tints[tint['id']] = tint
        else:
            re_dict = read_prog.match(line).groupdict()
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

def read_sequence(tints, read_files):
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
    name_to_idx=dict()
    for tidx,tint in enumerate(tints):
        for ridx,read in enumerate(tint['reads']):
            name_to_idx[read['name']] = (tidx,ridx)
            tints[tidx]['reads'][ridx]['seq']=''
    for read_file in read_files:
        read_file = open(read_file)
        first_line = read_file.readline()
        assert first_line[0] in ['@','>'], first_line
        if first_line[0] == '@':
            for idx,line in enumerate(chain([first_line], read_file)):
                if idx % 4 == 0:
                    assert line[0]=='@'
                    name = line[1:].split()[0]
                if idx % 4 == 1:
                    if name in name_to_idx:
                        tidx,ridx=name_to_idx[name]
                        tints[tidx]['reads'][ridx]['seq'] = line.rstrip()
                if idx % 4 == 2:
                    assert line.rstrip()=='+'
                if idx % 4 == 3:
                    if name in name_to_idx:
                        assert len(line.rstrip())==len(tints[tidx]['reads'][ridx]['seq']),line+'\n'+tints[tidx]['reads'][ridx]['seq']
        else:
            for line in chain([first_line], read_file):
                assert not line[0] in ['@','+']
                if line[0]=='>':
                    name = line[1:].split()[0]
                else:
                    if name in name_to_idx:
                        tidx,ridx=name_to_idx[name]
                        tints[tidx]['reads'][ridx]['seq'] += line.rstrip()
    for tint in tints:
        for read in tint['reads']:
            assert read['seq'] != ''
            read['length']=len(read['seq'])
            if read['strand']=='+':
                read['seq']=''.join(dna_id.get(x,'N') for x in read['seq'].upper())
            else:
                read['seq']=''.join(rev_comp.get(x,'N') for x in reversed(read['seq'].upper()))
    return

def forward_thread_cigar(cigar, t_goal, t_pos, q_pos):
    assert t_pos<=t_goal
    cig_idx = 0
    while t_pos < t_goal:
        c,t = cigar[cig_idx]
        # print(t_goal-t_pos, c,t)
        c = min(c, t_goal-t_pos)
        if t in ['M','X','=']:
            t_pos += c
            q_pos += c
        if t in ['D']:
            t_pos += c
        if t in ['I']:
            q_pos += c
        cig_idx+=1
    assert t_pos==t_goal
    return q_pos


def get_interval_start(start, read):
    '''
    Finds the first position in the read alignment intervals that aligns to start or after.
    Any (negative) offset is reported as slack.
    '''
    for (t_start,t_end,q_start,q_end,cigar) in read['intervals']:
        # We haven't passed the start location yet
        if t_end < start:
            continue
        # print((t_start, t_end),(q_start, q_end),'|',read['length'])
        if start < t_start:
            q_pos=q_start
            slack=start-t_start
        else:
            q_pos = forward_thread_cigar(cigar=cigar, t_goal=start, t_pos=t_start, q_pos=q_start)
            slack = 0
        assert slack<=0,(slack,t_start,start)
        assert q_start<=q_pos<=q_end, (q_start,q_pos,q_end)
        return q_pos,slack
    assert False

def get_interval_end(end, read):
    '''
    Finds the last location on the read that aligns on or before end.
    Any (negative) offset is reported as slack.
    '''
    for (t_start,t_end,q_start,q_end,cigar) in reversed(read['intervals']):
        # We haven't passed end yet
        if t_start > end:
            continue
        # print((t_start, t_end),(q_start, q_end),'|',read['length'])
        # the right most interval that covers end, ends before end
        if t_end<end:
            q_pos=q_end
            slack=t_end-end
        else:
            q_pos = forward_thread_cigar(cigar=cigar, t_goal=end, t_pos=t_start, q_pos=q_start)
            slack = 0
        assert slack<=0,(slack,t_end,end)
        assert 0<=q_pos<=q_end, (q_start,q_pos,q_end)
        return q_pos,slack
    assert False

def find_longest_poly(seq, match_score=1, mismatch_score=-2, char='A'):
    if len(seq)==0:
        return
    if seq[0]==char:
        scores=[match_score]
    else:
        scores=[0]
    for m in (match_score if c==char else mismatch_score for c in seq[1:]):
        scores.append(max(0,scores[-1]+m))
    for k,g in groupby(enumerate(scores),lambda x:x[1]>0):
        if not k:
            continue
        i,s = list(zip(*g))
        max_s,max_i=max(zip(s,i))
        l = max_i+1-i[0]
        yield i[0], l, seq[i[0]:i[0]+l].count(char)/l

def get_unaligned_gaps_and_polyA(read, segs):
    read['gaps']=set()
    if not 1 in read['data']:
        return
    intervals = list()
    for d,group in groupby(enumerate(read['data']), lambda x: x[1]):
        if d != 1:
            continue
        group = list(group)
        f_seg_idx = group[0][0]
        l_seg_idx = group[-1][0]
        intervals.append((f_seg_idx,l_seg_idx))
    assert len(intervals)>0, read['data']
    (f_seg_idx,_) = intervals[0]
    start = segs[f_seg_idx][0]
    q_ssc_pos,_ = get_interval_start(start=start, read=read)
    (_,l_seg_idx) = intervals[-1]
    end   = segs[l_seg_idx][1]
    q_esc_pos,_ = get_interval_end(end=end, read=read)
    assert 0<=q_ssc_pos<q_esc_pos<=read['length'], (q_ssc_pos,q_esc_pos,read['length'],start,end,read['intervals'],read['id'])
    s_polys = list()
    for char in ['A','T']:
        for i,l,p in find_longest_poly(read['seq'][0:q_ssc_pos], char=char):
            if l < 20 or p < 0.85:
                continue
            assert 0<=i<q_ssc_pos,(i,q_ssc_pos,read['length'])
            s_polys.append((i,l,p,char))
    if len(s_polys)>0:
        i,l,p,char = max(s_polys, key=lambda x: x[2])
        poly_to_gene_gap_size = q_ssc_pos-i-l
        assert 0<=poly_to_gene_gap_size<q_ssc_pos
        read['gaps'].add(
            'S{}_{}:{}'.format(char,l,poly_to_gene_gap_size)
        )
        read['gaps'].add(
            'SSC:{}'.format(i)
        )
    else:
        read['gaps'].add(
            'SSC:{}'.format(q_ssc_pos)
        )
    e_polys = list()
    for char in ['A','T']:
        for i,l,p in find_longest_poly(read['seq'][q_esc_pos:], char=char):
            if l < 20 or p < 0.85:
                continue
            assert 0<=i<read['length']-q_esc_pos,(i,q_esc_pos,read['length'])
            e_polys.append((i,l,p,char))
    if len(e_polys)>0:
        i,l,p,char = max(e_polys, key=lambda x: x[2])
        poly_to_gene_gap_size = i
        assert 0<=poly_to_gene_gap_size<read['length']-q_esc_pos, (q_esc_pos,i,l,p,read['length'],poly_to_gene_gap_size)
        read['gaps'].add(
            'E{}_{}:{}'.format(char,l,poly_to_gene_gap_size)
        )
        read['gaps'].add(
            'ESC:{}'.format(read['length']-q_esc_pos-poly_to_gene_gap_size)
        )
        assert read['length']-q_esc_pos-poly_to_gene_gap_size>0
    else:
        read['gaps'].add(
            'ESC:{}'.format(read['length']-q_esc_pos)
        )
    for i1,i2 in zip(intervals[:-1],intervals[1:]):
        (_,i1_l_seg_idx) = i1
        i1_end           = segs[i1_l_seg_idx][1]
        q_gap_start,start_slack      = get_interval_end(end=i1_end, read=read)
        (i2_f_seg_idx,_) = i2
        i2_start         = segs[i2_f_seg_idx][0]
        q_gap_end,end_slack        = get_interval_start(start=i2_start, read=read)
        assert 0<q_gap_start<=q_gap_end<read['length'],(q_gap_start,q_gap_end,read['length'])
        q_gap_size = q_gap_end-q_gap_start
        q_gap_size = max(0,q_gap_size+start_slack+end_slack)
        assert 0<=q_gap_size<read['length'],(q_gap_size,start_slack,end_slack)
        assert i1_l_seg_idx<i2_f_seg_idx
        read['gaps'].add(
            '{}-{}:{}'.format(i1_l_seg_idx,i2_f_seg_idx,q_gap_size),
        )
    read['gaps']=sorted(read['gaps'])

def optimize(candidate_y_idxs, C, start, end, smoothed_threshold, low, high, read_support):
    cov_mem = dict()
    yea_mem = dict()
    nay_mem = dict()
    amb_mem = dict()

    # print('Precomputing coverage mems for {}...'.format((start,end)))
    for i in range(start, end):
        for j in range(i, end+1):
            seg_len = (candidate_y_idxs[j]-candidate_y_idxs[i]+1)
            cov_mem[(i,j)] = (C[j]-C[i])/seg_len
            if seg_len < len(smoothed_threshold):
                h = smoothed_threshold[seg_len]
                l = 1-h
            else:
                h = high
                l = low
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

def run_optimize(candidate_y_idxs, fixed_c_idxs, coverage, smoothed_threshold, low_threshold, high_threshold, min_read_support_outside):
    final_c_idxs = set(fixed_c_idxs)
    for start,end in zip(fixed_c_idxs[:-1],fixed_c_idxs[1:]):
        D,B,max_d,max_b,in_mem,out_mem = optimize(
            candidate_y_idxs         = candidate_y_idxs,
            C                        = coverage,
            start                    = start,
            end                      = end,
            smoothed_threshold       = smoothed_threshold,
            low                      = low_threshold,
            high                     = high_threshold,
            read_support             = min_read_support_outside,
        )
        while max_b != (-1,-1,-1):
            final_c_idxs.update(max_b)
            # print('B:', max_b)
            # print('I:', in_mem[(max_b[0],max_b[1])])
            # print('O:', out_mem[max_b])
            max_b = B[max_b]
    return sorted(final_c_idxs)

def get_cumulative_coverage(candidate_y_idxs, y_idx_to_r_idxs):
    C = np.zeros((len(candidate_y_idxs)+1, y_idx_to_r_idxs.shape[1]), dtype=np.uint32)
    for C_idx,(cur_y_idx,nxt_y_idx) in enumerate(zip(candidate_y_idxs[:-1],candidate_y_idxs[1:]), start=1):
        C[C_idx] = y_idx_to_r_idxs[cur_y_idx:nxt_y_idx].sum(axis=0)
    for C_idx in range(1,len(C)):
        C[C_idx] += C[C_idx-1]
    return C

def segment(segment_args):
    tint, sigma, smoothed_threshold, high_threshold, variance_factor, max_candidates_per_seg, min_read_support_outside = segment_args
    low_threshold = 1 - high_threshold
    pos_to_Yy_idx = dict()
    Yy_idx_to_pos = list()
    Yy_idx_to_r_idxs = list()
    Y = list()
    # print('Generating index for tint {} with {} intervals'.format(tint['id'], tint['intervals']))
    for s,e in tint['intervals']:
        y_idx_to_pos = list()
        for p in range(s,e+1):
            assert not p in pos_to_Yy_idx
            pos_to_Yy_idx[p]=(len(Yy_idx_to_pos),len(y_idx_to_pos))
            y_idx_to_pos.append(p)
        Yy_idx_to_pos.append(y_idx_to_pos)
        Yy_idx_to_r_idxs.append(np.zeros((len(y_idx_to_pos), len(tint['reads'])), dtype=bool))
        Y.append(np.zeros(len(y_idx_to_pos)))
    assert len(pos_to_Yy_idx)==sum(len(y_idx_to_pos) for y_idx_to_pos in Y)
    # print('Building per read coverage data for tint {}'.format(tint['id']))
    for r_idx,read in enumerate(tint['reads']):
        read['data']=list()
        for ts,te,_,_,_ in read['intervals']:
            Y_idx,y_idx = pos_to_Yy_idx[ts]
            Y[Y_idx][y_idx] += 1
            Y_idx,y_idx = pos_to_Yy_idx[te]
            Y[Y_idx][y_idx] += 1
            for pos in range(ts,te):
                Y_idx,y_idx = pos_to_Yy_idx[pos]
                Yy_idx_to_r_idxs[Y_idx][y_idx][r_idx] = True

    Y = [gaussian_filter1d(y,sigma) for y in Y]
    Y_none_zero_vals = np.array([v for y in Y for v in y if v > 0])
    variance_threhsold = Y_none_zero_vals.mean() + variance_factor*Y_none_zero_vals.std()

    tint['final_positions'] = list()
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
        # print('Optimizing tint {} with:\n\tcandidate {} locations: {}\n\tfixed {} loations: {}'.format(tint['id'], len(candidate_y_idxs),candidate_y_idxs,len(fixed_c_idxs),[candidate_y_idxs[c_idx] for c_idx in fixed_c_idxs]))
        final_c_idxs=run_optimize(
            candidate_y_idxs=candidate_y_idxs,
            fixed_c_idxs=fixed_c_idxs,
            coverage=cumulative_coverage,
            smoothed_threshold=smoothed_threshold,
            low_threshold=low_threshold,
            high_threshold=high_threshold,
            min_read_support_outside=min_read_support_outside
        )
        final_y_idxs = [candidate_y_idxs[c_idx] for c_idx in final_c_idxs]
        tint['final_positions'].extend([Yy_idx_to_pos[Y_idx][y_idx] for y_idx in final_y_idxs])
        for r_idx,read in enumerate(tint['reads']):
            for s,e in zip(final_c_idxs[:-1],final_c_idxs[1:]):
                seg_len = e-s+1
                if seg_len < len(smoothed_threshold):
                    ht = high_threshold
                    lt = low_threshold
                else:
                    ht = smoothed_threshold[seg_len]
                    lt = 1-ht
                cov_ratio = (cumulative_coverage[e][r_idx]-cumulative_coverage[s][r_idx])/seg_len
                if cov_ratio > ht:
                    read['data'].append(1)
                elif cov_ratio < lt:
                    read['data'].append(0)
                else:
                    read['data'].append(2)
            read['data'].append(0)
    tint['segs']=[(s,e) for s,e in zip(tint['final_positions'][:-1],tint['final_positions'][1:])]
    # print('Extracting unaligned gaps and polyA tail data from reads for tint {}'.format(tint['id']))
    for read in tint['reads']:
        read['data'].pop()
        assert len(read['data'])==len(tint['segs']),(read['data'],tint['segs'])
        get_unaligned_gaps_and_polyA(read=read, segs=tint['segs'])
    return tint

def smooth_threshold(threshold):
    smooth = list()
    while True:
        x = len(smooth)
        y = threshold/(1 + ((threshold-.5)/.5)*exp(-0.05*x))
        if x>5 and x*(threshold-y)<0.5:
            break
        smooth.append(round(y,2))
        assert len(smooth)<1000
    return smooth

def main():
    args = parse_args()

    tints = read_split(args.split_tsv)
    read_sequence(tints, args.reads)

    segment_args = list()
    for tint in tints[:]:
        segment_args.append((
            tint,
            args.sigma,
            smooth_threshold(threshold=args.threshold),
            args.threshold,
            args.variance_factor,
            args.max_candidates_per_seg,
            args.min_read_support_outside,
        ))
    out_file = open(args.output, 'w')
    with Pool(args.threads) as p:
        if args.threads == 1:
            p.close()
        for idx,tint in enumerate(p.imap_unordered(segment, segment_args, chunksize=20)) if args.threads>1 else enumerate(map(segment, segment_args)):
            tints[tint['id']]=tint
            record = list()
            record.append('#{}'.format(tint['chr']))
            record.append(str(tint['id']))
            record.append(','.join(map(str,tint['final_positions'])))
            print('\t'.join(record), file=out_file)
            for read in tint['reads']:
                record = list()
                record.append(str(read['id']))
                record.append(read['name'])
                record.append(read['chr'])
                record.append(read['strand'])
                record.append(str(read['tint']))
                record.append(''.join(map(str,read['data'])))
                record.append(''.join('{},'.format(g) for g in read['gaps']))
                print('\t'.join(record), file=out_file)
            print('Done with {}-th transcriptional multi-intervals ({}/{})'.format(tint['id'], idx+1,len(tints)))
    out_file.close()

if __name__ == "__main__":
    main()
