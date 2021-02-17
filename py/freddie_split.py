#!/usr/bin/env python3
import argparse
import os
import resource
import re
from operator import itemgetter
from collections import deque
from multiprocessing import Pool

import pysam

cigar_re = re.compile(r'(\d+)([M|I|D|N|S|H|P|=|X]{1})')

query_consuming = [
    pysam.CINS,
    pysam.CSOFT_CLIP,
    pysam.CMATCH,
    pysam.CEQUAL,
    pysam.CDIFF,
]
target_consuming = [
    pysam.CDEL,
    pysam.CMATCH,
    pysam.CEQUAL,
    pysam.CDIFF,
]
exon_consuming = [
    pysam.CINS,
    pysam.CDEL,
    pysam.CMATCH,
    pysam.CEQUAL,
    pysam.CDIFF,
]
intron_consuming = [
    pysam.CINS,
    pysam.CDEL,
    pysam.CMATCH,
    pysam.CEQUAL,
    pysam.CDIFF,
]
target_and_query_consuming = [
    pysam.CMATCH,
    pysam.CEQUAL,
    pysam.CDIFF,
]
target_skipping = [
    pysam.CDEL,
    pysam.CREF_SKIP,
]
cop_to_str = [
    'M',
    'I',
    'D',
    'N',
    'S',
    'H',
    'P',
    '=',
    'X',
    'B',
]


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
    parser.add_argument("-t",
                        "--threads",
                        default=1,
                        type=int,
                        help="Number of threads. Max # of threads used is # of contigs. Default: 1")
    parser.add_argument("-o",
                        "--outdir",
                        type=str,
                        default='freddie_split/',
                        help="Path to output directory. Default: freddie_split/")
    args = parser.parse_args()
    assert args.threads > 0
    return args


# def fix_cigar(cigar):
#     fixed_cigar = list()
#     fixed_cigar.append(cigar[0])
#     for t, c in cigar[1:]:
#         last_t, last_c = fixed_cigar[-1]
#         if t == last_t:
#             fixed_cigar[-1] = (t, last_c+c)
#             continue
#         if t in target_skipping and last_t in target_skipping:
#             fixed_cigar[-1] = (pysam.CREF_SKIP, c+last_c)
#             continue
#         fixed_cigar.append((t, c))
#     return fixed_cigar


def fix_intervals(intervals):
    for ts, te, qs, qe, cigar in intervals:
        if len(cigar) == 0:
            continue
        (t, c) = cigar[0]
        if t == pysam.CDEL:
            ts += c
            cigar = cigar[1:]
        if len(cigar) == 0:
            continue
        (t, c) = cigar[-1]
        if t == pysam.CDEL:
            te -= c
            cigar = cigar[:-1]
        if ts < te:
            yield (ts, te, qs, qe, cigar)


def get_intervals(aln):
    cigar = aln.cigartuples
    qstart = 0
    if cigar[0][0] == pysam.CSOFT_CLIP:
        qstart += cigar[0][1]
    qlen = 0
    for t, c in cigar:
        if t in query_consuming:
            qlen += c
    assert qlen == len(aln.query_sequence)
    qend = qlen
    if cigar[-1][0] == pysam.CSOFT_CLIP:
        qend -= cigar[-1][1]
    assert qend > qstart
    tstart = aln.reference_start
    tend = tstart + 1

    for t, c in cigar:
        assert(0 <= t < 10)
    for t, c in cigar:
        if t in target_consuming:
            tend += c
    qstart_c = qstart
    qend_c = qstart
    tstart_c = tstart
    tend_c = tstart

    intervals = list()
    interval_cigar = list()
    for t, c in cigar:
        if t in exon_consuming:
            interval_cigar.append((t, c))
        if t == pysam.CDEL:
            tend_c += c
        if t == pysam.CINS:
            qend_c += c
        if t in target_and_query_consuming:
            tend_c += c
            qend_c += c
        if t == pysam.CREF_SKIP:
            intervals.append((
                tstart_c,
                tend_c,
                qstart_c,
                qend_c,
                interval_cigar,
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
            interval_cigar,
        ))
    return list(fix_intervals(intervals))


def read_sam(sam, contig):
    reads = list()
    start, end = None, None
    for aln in sam.fetch(contig=contig):
        if aln.is_unmapped or aln.is_supplementary or aln.is_secondary or aln.reference_name == None:
            continue
        assert aln.reference_name == contig, '{} : {}'.format(
            aln.reference_name, contig)
        read = dict(
            id=len(reads),
            name=aln.query_name,
            contig=aln.reference_name,
            strand='-' if aln.is_reverse else '+',
            simple_tints=list(),
            tint=None,
            intervals=[(st, et, sr, er, c) for (st, et, sr, er, c)
                       in get_intervals(aln) if st != et and sr != er],
        )
        s, _, _, _, _ = read['intervals'][0]
        _, e, _, _, _ = read['intervals'][-1]
        if (start, end) == (None, None):
            start, end = s, e
        if s > end:
            yield reads
            reads = list()
            read['id'] = len(reads)
            end = e
        end = max(end, e)
        reads.append(read)
    if len(reads) > 0:
        yield reads


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
    print('[freddie_split] Splitting reads...')
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
        print('[freddie_split.py] Sorting contig {}...'.format(c))
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


def run_split(split_args):
    bam, contig, outdir = split_args
    sam = pysam.AlignmentFile(bam, 'rb')
    rname_to_tint = dict()
    tint_id = 0
    contig_outdir = '{}/{}'.format(outdir, contig)
    for reads in read_sam(sam=sam, contig=contig):
        tints = get_transcriptional_intervals(reads=reads)
        for tint in tints:
            if tint_id == 0:
                print('[freddie_split.py] Splitting contig {}'.format(contig))
                os.makedirs(contig_outdir, exist_ok=False)
            write_tint(contig_outdir, contig, tint_id,
                       tint, reads, rname_to_tint)
            tint_id += 1
    return contig, rname_to_tint


def write_tint(contig_outdir, contig, tint_id, tint, reads, rname_to_tint):
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
            record.append(parse_interval_field(interval))
        outfile.write('\t'.join(record))
        outfile.write('\n')
    outfile.close()


def parse_interval_field(interval):
    return '{}-{}:{}-{}:{}'.format(interval[0], interval[1], interval[2], interval[3], ''.join(('{}{}'.format(c, cop_to_str[t]) for t, c in interval[4])))


def main():
    args = parse_args()

    args.outdir = args.outdir.rstrip('/')
    os.makedirs(args.outdir, exist_ok=True)
    print('[freddie_split.py] Running split with args:', args)

    contigs = {x['SN']: x['LN']
               for x in pysam.AlignmentFile(args.bam, 'rb').header['SQ']}
    args.threads = min(args.threads, len(contigs))
    split_args = list()
    for contig, _ in sorted(contigs.items(), key=itemgetter(1), reverse=True):
        split_args.append((
            args.bam,
            contig,
            args.outdir,
        ))
    rname_to_tint = dict()
    final_contigs = list()
    if args.threads > 1:
        p = Pool(args.threads)
        for contig, rname_to_tint_thread in p.imap_unordered(run_split, split_args, chunksize=1):
            if len(rname_to_tint_thread) == 0:
                continue
            print('[freddie_split] Done with contig {}'.format(contig))
            rname_to_tint = {**rname_to_tint, **rname_to_tint_thread}
            final_contigs.append(contig)
        p.close()
    else:
        for contig, rname_to_tint_thread in map(run_split, split_args):
            if len(rname_to_tint_thread) == 0:
                continue
            print('[freddie_split] Done with contig {}'.format(contig))
            rname_to_tint = {**rname_to_tint, **rname_to_tint_thread}
            final_contigs.append(contig)

    RLIMIT_NOFILE_soft,RLIMIT_NOFILE_hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    if RLIMIT_NOFILE_hard < len(contigs)+10:
        assert False, 'Number of contigs in reference is much larger than system hard limit on open files! {} vs {}'.format(len(contigs), RLIMIT_NOFILE_hard)
    if RLIMIT_NOFILE_soft < len(contigs)+10:
        resource.setrlimit(resource.RLIMIT_NOFILE, (len(contigs)+10, RLIMIT_NOFILE_hard))
    split_reads(read_files=args.reads,
                rname_to_tint=rname_to_tint, contigs=final_contigs, outdir=args.outdir)

if __name__ == "__main__":
    main()
