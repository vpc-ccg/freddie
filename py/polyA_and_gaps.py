#!/usr/bin/env python3
import argparse

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

    segs = [int(x.rstrip()) for x in open(args.segs)]
    segs = [(s,e) for s,e in zip(segs[:-1],segs[1:])]
    for s,e in segs:
        assert s<e
    pos_to_rid,rid_to_len,rid_to_intervals,read_name_to_id,t_len = read_paf(args.paf)
    for name,rid in {name.rstrip():rid for (rid,name) in enumerate(open(args.names))}.items():
        assert read_name_to_id[name]==rid
    rid_to_data = [[x for x in d.rstrip()] for d in open(args.data)]
    assert len(rid_to_data)==len(read_name_to_id)==len(rid_to_len)==len(rid_to_intervals)
    for d in rid_to_data:
        assert len(d)==len(segs)

    # out_file = open('{}.gaps'.format(args.out_prefix), 'w+')
    # rid_to_unaln_gaps = get_unaligned_gaps(data=data,brks=peak_positions,t_len=t_len,rid_to_intervals=rid_to_intervals,rid_to_len=rid_to_len)
    # for rid in range(len(data)):
    #     print('\t'.join(['{}-{}-{}'.format(*x) for x in rid_to_unaln_gaps[rid]]),file=out_file)
    # out_file.close()

if __name__ == "__main__":
    main()
