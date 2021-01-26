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
    parser.add_argument("-b",
                        "--bam",
                        type=str,
                        required=True,
                        help="Path to sorted and indexed BAM file of reads. Assumes splice aligner is used to the genome. Prefers deSALT")
    parser.add_argument("-r",
                        "--reads",
                        nargs="+",
                        type=str,
                        required=True,
                        help="Space separated paths to reads in FASTQ or FASTA format used to extract polyA tail information")
    parser.add_argument("-o",
                        "--outdir",
                        type=str,
                        default='freddie_split/',
                        help="Path to output directory. Default: freddie_split/")
    args = parser.parse_args()
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
    start, end = None, None
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
        s, e, _, _, _ = reads[-1]['intervals'][-1]
        if (start, end) == (None, None):
            start, end = s, e
        if s > end:
            yield reads
            reads = list()
            start, end = s, e
        if e > end:
            end = e


def get_transcriptional_intervals(reads):
    intervals = list()
    start, end = None, None
    rids = list()
    for s, e, rid in sorted((i[0], i[1], read['id']) for read in reads for i in read['intervals']):
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
        if len(rids) < 3:
            continue
        for rid in rids:
            reads[rid]['tint'] = len(multi_tints)
        multi_tints.append(dict(intervals=sorted(
            group_intervals), rids=sorted(rids)))
    assert all(enqueued)
    return multi_tints


def split_reads(read_files, rname_to_tint, contigs, outdir):
    outfiles = {c: open('{}/{}/reads.tsv'.format(outdir, c), 'w+')
                for c in contigs}
    for read_file in read_files:
        for idx, line in enumerate(open(read_file)):
            if idx == 0:
                if line[0] == '@':
                    mod = 4
                elif line[0] == '>':
                    mod = 2
                else:
                    assert False, 'Invalid fasta/q file ' + read_file
            if idx % mod == 0:
                rname = line.rstrip().split()[0][1:]
                (contig, tint_id, rid) = rname_to_tint.get(
                    rname, (None, None, None))
            if idx % mod == 1 and contig != None:
                seq = line.rstrip()
                record = list()
                record.append(str(rid))
                record.append(str(contig))
                record.append(str(tint_id))
                record.append(seq)
                outfiles[contig].write('\t'.join(record))
                outfiles[contig].write('\n')
    for outfile in outfiles.values():
        outfile.close()
    for c in contigs:
        path = '{od}/{c}/reads.tsv'.format(od=outdir, c=c)
        os.system(
            'sort -k3,3n {} > {}_sorted'.format(path, path))
        os.system('mv {}_sorted {}'.format(path, path))
        last_tint = None
        for line in open(path):
            rid, contig, tint_id, _ = line.rstrip().split('\t')
            if last_tint == None:
                last_tint = tint_id
                outfile = open(
                    '{}/{}/reads_{}_{}.tsv'.format(outdir, c, c, tint_id), 'w+')
            if last_tint != tint_id:
                outfile.close()
                last_tint = tint_id
                outfile = open(
                    '{}/{}/reads_{}_{}.tsv'.format(outdir, c, c, tint_id), 'w+')
            outfile.write(line)
        outfile.close()
        # os.remove(path)


def run_split(sam, reads, outdir):
    sam = pysam.AlignmentFile(sam, 'rb')

    rname_to_tint = dict()
    contigs = {x['SN'] for x in sam.header['SQ']}
    tint_id = 0
    for contig in contigs:
        print('[freddie_split.py] Splitting contig {}'.format(contig))
        for reads in read_sam(sam=sam, contig=contig):
            contig_outdir = '{}/{}'.format(outdir, contig)
            os.makedirs(contig_outdir, exist_ok=False)
            tints = get_transcriptional_intervals(reads=reads)
            for tint in tints:
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
                    rname_to_tint[read['name']] = (contig, tint_id, rid)
                    record = list()
                    record.append(str(read['id']))
                    record.append(read['name'])
                    record.append(read['contig'])
                    record.append(read['strand'])
                    record.append(str(tint_id))
                    for interval in read['intervals']:
                        record.append('{}-{}:{}-{}:{}'.format(*interval))
                    outfile.write('\t'.join(record))
                    outfile.write('\n')
                outfile.close()
                tint_id += 1
    split_reads(read_files=reads,
                rname_to_tint=rname_to_tint, contigs=contigs, outdir=outdir)


def main():
    args = parse_args()

    args.outdir = args.outdir.rstrip('/')
    os.makedirs(args.outdir, exist_ok=True)
    print('[freddie_split.py] Running split with args:', args)
    run_split(sam=args.sam, reads=args.reads, outdir=args.outdir)


if __name__ == "__main__":
    main()
