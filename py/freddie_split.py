#!/usr/bin/env python3
import argparse
import os
import re
import pysam
from collections import deque

cigar_re = re.compile(r'(\d+)([M|I|D|N|S|H|P|=|X]{1})')


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract alignment information from BAM/SAM file and splits reads into distinct transcriptional intervals")
    parser.add_argument("-s",
                        "--sam",
                        type=str,
                        required=True,
                        help="Path to BAM/SAM file of reads. Assumes splice aligner is used to the genome. Prefers deSALT")
    parser.add_argument("-f",
                        "--sam-format",
                        type=str,
                        default=None,
                        help="SAM file format: bam or sam. Default: infers format from --sam file extension")
    parser.add_argument("-o",
                        "--outdir",
                        type=str,
                        default='freddie_split/',
                        help="Path to output directory. Default: freddie_split/")
    args = parser.parse_args()
    assert args.sam_format in ['sam', 'bam', None]
    return args


def fix_cigar(cigar):
    fixed_cigar = list()
    fixed_cigar.append(cigar[0])
    for c, t in cigar[1:]:
        last_c, last_t = fixed_cigar[-1]
        if t == last_t:
            fixed_cigar[-1] = (last_c+c, t)
            continue
        if t in ['D', 'N'] and last_t in ['D', 'N']:
            fixed_cigar[-1] = (c+last_c, 'N')
            continue
        fixed_cigar.append((c, t))
    return fixed_cigar


def get_intervals(aln):
    if aln.cigarstring == None:
        return list()
    cigar = [x for x in cigar_re.findall(aln.cigarstring)]
    assert sum(len(x[0])+len(x[1]) for x in cigar) == len(
        aln.cigarstring), 'Errorenous cigar for aln: {}'.format(aln)
    cigar = [(int(x[0]), x[1]) for x in cigar_re.findall(aln.cigarstring)]
    qstart = 0
    if cigar[0][1] in ['S']:
        qstart += cigar[0][0]
    qlen = 0
    for c, t in cigar:
        if t in ['I', 'S', 'M', '=', 'X']:
            qlen += c
    assert qlen == len(aln.query_sequence)
    qend = qlen
    if cigar[-1][1] in ['S']:
        qend -= cigar[-1][0]
    assert qend > qstart
    tstart = aln.reference_start
    tend = tstart + 1

    cigar = fix_cigar(cigar)
    for c, t in cigar:
        if t in ['D', 'M', '=', 'X']:
            tend += c
    qstart_c = qstart
    qend_c = qstart
    tstart_c = tstart
    tend_c = tstart

    intervals = list()
    interval_cigar = list()
    for c, t in cigar:
        if t in ['D', 'I', 'X', '=', 'M']:
            interval_cigar.append('{}{}'.format(c, t))
        if t in ['D']:
            tend_c += c
        if t in ['I']:
            qend_c += c
        if t in ['X', '=', 'M']:
            tend_c += c
            qend_c += c
        if t == 'N':
            intervals.append((
                tstart_c,
                tend_c,
                qstart_c,
                qend_c,
                ''.join(interval_cigar)
            ))
            interval_cigar = list()
            tend_c += c
            tstart_c = tend_c
            qstart_c = qend_c
    if tstart_c < tend_c:
        intervals.append((
            tstart_c,
            tend_c,
            qstart_c,
            qend_c,
            ''.join(interval_cigar)
        ))
    return intervals


def read_sam(sam, contig):
    reads = list()
    for aln in sam.fetch(contig=contig):
        if aln.is_supplementary or aln.is_secondary or aln.reference_name == None:
            continue
        assert aln.reference_name == contig, '{} : {}'.format(
            aln.reference_name, contig)
        reads.append(dict(
            id=len(reads),
            name=aln.query_name,
            contig=aln.reference_name,
            strand='-' if aln.is_reverse else '+',
            simple_tints=list(),
            tint=None,
            intervals=[(st, et, sr, er, c) for (st, et, sr, er, c)
                       in get_intervals(aln) if st != et and sr != er],
        ))
    return reads


def get_transcriptional_intervals(reads, contig):
    intervals = list()
    start, end = None, None
    rids = list()
    for s, e, rid in sorted((i[0], i[1], read['id']) for read in reads if read['contig'] == contig for i in read['intervals']):
        if (start, end) == (None, None):
            start, end = s, e
        if s > end:
            intervals.append(dict(
                tint=-1,
                start=start,
                end=end,
                rids=rids,
            ))
            start = s
            end = e
            rids = list()
        assert start <= s
        end = max(end, e)
        rids.append(rid)
        reads[rid]['simple_tints'].append(len(intervals))
    if (start, end) == (None, None):
        return list()
    intervals.append(dict(
        tint=-1,
        start=start,
        end=end,
        rids=rids,
    ))

    enqueued = [False for _ in intervals]
    multi_tints = list()
    for idx in range(len(intervals)):
        if enqueued[idx]:
            continue
        group = list()
        queue = deque()
        queue.append(idx)
        enqueued[idx] = True
        while len(queue) > 0:
            tint = queue.pop()
            group.append(tint)
            for rid in intervals[tint]['rids']:
                for i in reads[rid]['simple_tints']:
                    if not enqueued[i]:
                        enqueued[i] = True
                        queue.append(i)
        rids = set()
        group_intervals = list()
        for tint in group:
            rids.update(intervals[tint]['rids'])
            group_intervals.append(
                (intervals[tint]['start'], intervals[tint]['end']))
        for rid in rids:
            reads[rid]['tint'] = len(multi_tints)
        multi_tints.append(dict(intervals=sorted(
            group_intervals), rids=sorted(rids)))
    assert all(enqueued)
    return multi_tints


def run_split(args):

    if args.sam_format == 'sam':
        sam = pysam.AlignmentFile(sam, 'r')
    elif args.sam_format == 'bam':
        sam = pysam.AlignmentFile(sam, 'rb')

    contigs = {x['SN']: int(x['LN']) for x in sam.header['SQ']}

    for contig in contigs.keys():
        print('[freddie_split.py] Splitting contig {}', contig)
        reads = read_sam(sam=sam, contig=contig)

        contig_outdir = '{}/{}'.format(args.outdir, contig)
        os.makedirs(contig_outdir, exist_ok=False)
        tints = get_transcriptional_intervals(reads=reads, contig=contig)
        for tint_id, tint in enumerate(tints):
            if (100*tint_id)//len(tints) % 10 == 0:
                print()
                print('[freddie_split.py] Done with {}/{}', tint_id, len(tints))
            outfile = open(
                '{}/split_{}_{}.tsv'.format(contig_outdir, contig, tint_id), 'w+')

            record = list()
            record.append('#{}'.format(contig))
            record.append('{}'.format(tint_id))
            record.append(','.join('{}-{}'.format(s, e)
                                   for s, e in tint['intervals']))
            record.append(str(len(tint['rids'])))
            outfile.write('\t'.join(record))
            outfile.write('\n')
            for rid in tint['rids']:
                read = reads[rid]
                record = list()
                record.append(str(read['id']))
                record.append(read['name'])
                record.append(read['contig'])
                record.append(read['strand'])
                record.append(tint_id)
                for interval in read['intervals']:
                    record.append('{}-{}:{}-{}:{}'.format(*interval))
                outfile.write('\t'.join(record))
                outfile.write('\n')
            outfile.close()


def main():
    args = parse_args()

    if args.sam_format == None:
        args.sam_format = 'sam'
        if args.sam.endswith('.bam'):
            args.sam_format = 'bam'
    args.outdir = args.outdir.rstrip('/')
    os.makedirs(args.outdir, exist_ok=True)
    print('[freddie_split.py] Running split with args:', args)
    run_split(args=args)


if __name__ == "__main__":
    main()
