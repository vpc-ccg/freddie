#!/usr/bin/env python3
import argparse
import re
from itertools import groupby

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
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to PAF file of read alignments")
    parser.add_argument("-q",
                        "--fastq",
                        type=str,
                        required=True,
                        help="Path to FASTQ file of read alignments")
    parser.add_argument("-s",
                        "--segs",
                        type=str,
                        required=True,
                        help="Path to segment breakpoints TXT file")
    parser.add_argument("-d",
                        "--data",
                        type=str,
                        required=True,
                        help="Path to DATA file")
    parser.add_argument("-n",
                        "--names",
                        type=str,
                        required=True,
                        help="Path to names TXT file")
    parser.add_argument("-op",
                        "--out-prefix",
                        type=str,
                        required=True,
                        help="Output prefix that does not include .TXT part")
    args = parser.parse_args()
    return args

def read_paf(paf, range_len=0):
    is_first = True
    reads = list()
    read_name_to_id = dict()
    for line in open(paf):
        line = line.rstrip().split('\t')
        if is_first:
            t_len = int(line[6])
            t_name = line[5]
            is_first = False
        assert t_len == int(line[6]) and t_name == line[5], 'Multiple targets detected in PAF file!\n{}'.format(line)
        if not any('oc:c:1' in tag for tag in line[12:]):
            continue
        name = line[0]
        if not name in read_name_to_id:
            read_name_to_id[name] = len(reads)
            reads.append(dict(
                name=name,
                intervals=list(),
                length=int(line[1]),
                strand=line[4],
            ))
        rid = read_name_to_id[name]
        assert reads[rid]['length'] == int(line[1])
        assert reads[rid]['strand'] == line[4]
        cigar = None
        for tag in line[12:]:
            if tag[:len('cg:Z:')]=='cg:Z:':
                assert cigar == None
                cigar = [(int(x[0]),x[1]) for x in re.findall(r'(\d+)([M|I|D|N|S|H|P|=|X]{1})', tag[len('cg:Z:'):])]
                assert sum([len(str(x[0])+x[1]) for x in cigar])==len(tag[len('cg:Z:'):]),'Something wrong with line:\n{}'.format(line)
        assert not cigar == None, 'PAF record has no cigar:\n{}'.format(line)
        t_start = int(line[7])
        t_end   = int(line[8])
        q_start = int(line[2])
        q_end   = int(line[3])
        assert 0 <= t_start < t_end <= t_len
        assert 0 <= q_start < q_end <= reads[rid]['length'], '{}<{}<{}:\n{}'.format(q_start,q_end,reads[rid]['length'],line)
        if reads[rid]['strand']=='-':
            q_end,q_start = reads[rid]['length']-q_start-1,reads[rid]['length']-q_end-1
        t_interval = (t_start, t_end)
        q_interval = (q_start, q_end)
        reads[rid]['intervals'].append(((t_interval, q_interval), cigar))
    for read in reads:
        read['intervals'].sort()
        for ((ti_1,qi_1),_),((ti_2,qi_2),_) in zip(read['intervals'][:-1],read['intervals'][1:]):
            assert ti_1[1]<=ti_2[0]
            assert qi_1[1]<=qi_2[0]
    return reads,read_name_to_id,t_len

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
    return q_pos


def get_interval_start(start, end, read):
    for ((t_start, t_end),(q_start, q_end)),cigar in read['intervals']:
        if t_end < start or t_start>end:
            continue
        # print((t_start, t_end),(q_start, q_end),'|',read['length'])
        if start < t_start:
            q_pos=q_start
            slack=t_start-start
        else:
            q_pos = forward_thread_cigar(cigar=cigar, t_goal=start, t_pos=t_start, q_pos=q_start)
            slack = 0
        assert q_pos<q_end, (q_pos,q_end)
        return q_pos
    assert False

def get_interval_end(start, end, read):
    for ((t_start, t_end),(q_start, q_end)),cigar in reversed(read['intervals']):
        if t_start > end or t_end<start:
            continue
        # print((t_start, t_end),(q_start, q_end),'|',read['length'])
        if t_end<end:
            q_pos=q_start
            slack=end-t_end
        else:
            q_pos = forward_thread_cigar(cigar=cigar, t_goal=end-1, t_pos=t_start, q_pos=q_start)
            slack = 0
        assert q_pos<q_end, (q_pos,q_end)
        return q_pos
    assert False

def get_unaligned_gaps(reads, segs, tlen):
    for read in reads:
        read['gaps']=set()
        if not 1 in read['data']:
            continue
        intervals = list()
        for d,group in groupby(enumerate(read['data']), lambda x: x[1]):
            if d != 1:
                continue
            group = list(group)
            f_seg_idx = group[0][0]
            l_seg_idx = group[-1][0]
            intervals.append((f_seg_idx,l_seg_idx))
        assert len(intervals)>0, read['data']

        for i1,i2 in zip(intervals[:-1],intervals[1:]):
            (i1_f_seg_idx,i1_l_seg_idx) = i1
            i1_start = segs[i1_f_seg_idx][0]
            i1_end   = segs[i1_l_seg_idx][1]
            q_gap_start = get_interval_end(start=i1_start, end=i1_end, read=read)
            (i2_f_seg_idx,i2_l_seg_idx) = i2
            i2_start = segs[i2_f_seg_idx][0]
            i2_end   = segs[i2_l_seg_idx][1]
            q_gap_end = get_interval_start(start=i2_start, end=i2_end, read=read)
            q_gap_size = q_gap_end-q_gap_start
            assert q_gap_size>0
            assert i1_l_seg_idx<i2_f_seg_idx
            read['gaps'].add((i1_l_seg_idx,i2_f_seg_idx,q_gap_size))
def get_read_seqs(reads, read_name_to_id, fastq):
    for idx,line in enumerate(open(fastq)):
        if idx%4==0:
            rid = read_name_to_id[line.rstrip().split()[0][1:]]
        if idx%4==1:
            reads[rid]['seq'] = line.rstrip()
    for read in reads:
        seq=read['seq']
        assert len(seq)==read['length'], (read['name'],len(seq),seq,read['length'])
        if read['strand']=='+':
            read['seq']=''.join(dna_id.get(x,'N') for x in seq.upper())
        else:
            read['seq']=''.join(rev_comp.get(x,'N') for x in reversed(seq.upper()))

def main():
    args = parse_args()

    segs = [int(x.rstrip()) for x in open(args.segs)]
    segs = [(s,e) for s,e in zip(segs[:-1],segs[1:])]
    for s,e in segs:
        assert s<e
    reads,read_name_to_id,tlen = read_paf(args.paf)
    for name,rid in {name.rstrip():rid for (rid,name) in enumerate(open(args.names))}.items():
        assert read_name_to_id[name]==rid
        assert reads[rid]['name']==name
    for rid,d in enumerate(open(args.data)):
        reads[rid]['data']=[int(x) for x in d.rstrip()]
        assert len(reads[rid]['data'])==len(segs)
    get_read_seqs(reads=reads, read_name_to_id=read_name_to_id, fastq=args.fastq)
    get_unaligned_gaps(reads=reads, segs=segs, tlen=tlen)

    # out_file = open('{}.gaps'.format(args.out_prefix), 'w+')
    # rid_to_unaln_gaps = get_unaligned_gaps(data=data,brks=peak_positions,t_len=t_len,rid_to_intervals=rid_to_intervals,rid_to_len=rid_to_len)
    # for rid in range(len(data)):
    #     print('\t'.join(['{}-{}-{}'.format(*x) for x in rid_to_unaln_gaps[rid]]),file=out_file)
    # out_file.close()

if __name__ == "__main__":
    main()
