#!/usr/bin/env python3
import argparse
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(
        description="Output coverage per simulated transcript")
    parser.add_argument("-g",
                        "--annotation-gtf",
                        type=str,
                        required=True,
                        help="Path to annotation GTF file of the transcripts")
    parser.add_argument("-s",
                        "--split-tsv",
                        type=str,
                        required=True,
                        help="Path to Freddie split TSV file of the reads")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default='freddie_cov.tsv',
                        help="Path to output file. Default: freddie_cov.tsv")
    parser.add_argument("-prt",
                        "--per-read-threshold",
                        type=float,
                        default=0.0,
                        help="Percentage of transcript length required to count a read simualated from that transcript as covering it. Default: 0.0")
    parser.add_argument("-ptt",
                        "--per-transcript-threshold",
                        type=int,
                        default=3,
                        help="Minimum reads coverage for a transcript to output. Default: 3")
    args = parser.parse_args()
    assert(args.per_transcript_threshold >= 0)
    assert(1 >= args.per_read_threshold >= 0)
    return args

def read_gtf(annotation_gtf):
    tids = {'others':dict(chrom='None', intervals=list(), reads=list())}
    for line in open(annotation_gtf):
        if line[0]=='#':
            continue
        line=line.rstrip().split('\t')
        if line[2] == 'transcript':
            tid = [x.strip(' "').split('"')[1] for x in line[8].strip(';').split(';') if x.startswith(' transcript_id')][0]
            tids[tid] = dict(
                chrom=line[0],
                intervals=list(),
                reads=list(),
            )
        if line[2] == 'exon':
            assert 'transcript_id "{}";'.format(tid) in line[8]
            tids[tid]['intervals'].append((
                int(line[3])-1,
                int(line[4]),
            ))
    return tids

def read_split(split_tsv, tids):
    for line in open(split_tsv):
        if line[0]=='#':
            continue
        line=line.rstrip().split('\t')
        tid=line[1].split('_')[0]
        if not tid in tids:
            tid='others'
        tids[tid]['reads'].append(dict(
            chrom=line[2],
            intervals=[tuple(map(int,x.split(':')[0].split('-'))) for x in line[5:]]
        ))

def overlap(chrom, intervals, queries):
    if len(intervals) == 0:
        return [0.0 for _ in queries]
    f=min(i[0] for i in intervals)
    l=max(i[1] for i in intervals)
    A = np.array([False]*(l-f), dtype=bool)
    for s,e in intervals:
        A[s-f:e-f]=True
    d = sum(A)
    result=list()
    for q in queries:
        o = 0.0
        if q['chrom']==chrom:
            o = sum(sum(A[s-f:e-f]) for s,e in q['intervals'])
        result.append(o/d)
    return result

def output_coverage(tids, per_read_threshold, per_transcript_threshold, outpath):
    out_file = open(outpath, 'w+')
    for tid,v in tids.items():
        s = sum(o >= per_read_threshold for o in overlap(v['chrom'], v['intervals'], v['reads']))
        if s >= per_transcript_threshold:
            print('{}\t{}\t'.format(tid,len(v['reads']), s), file=out_file)
    out_file.close()

def main():
    args = parse_args()

    tids = read_gtf(args.annotation_gtf)
    read_split(args.split_tsv, tids)
    output_coverage(tids, args.per_read_threshold, args.per_transcript_threshold, args.output)

if __name__ == "__main__":
    main()
