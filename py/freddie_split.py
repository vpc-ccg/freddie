#!/usr/bin/env python3
import argparse
import os
import resource
import re
from operator import itemgetter
from collections import deque
from multiprocessing import Pool
import gzip

import pysam
import networkx as nx
import numpy as np

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
                        help="Space separated paths to reads in FASTQ or FASTA format used to extract polyA tail information. If the file ends with .gz, it will be read using gzip")
    parser.add_argument("-t",
                        "--threads",
                        default=1,
                        type=int,
                        help="Number of threads. Max # of threads used is # of contigs. Default: 1")
    parser.add_argument("--contig-min-size",
                        default=1_000_000,
                        type=int,
                        help="Minimum contig size. Any contig with less size will not be processes. Default: 1,000,000")
    parser.add_argument("-o",
                        "--outdir",
                        type=str,
                        default='freddie_split/',
                        help="Path to output directory. Default: freddie_split/")
    args = parser.parse_args()
    assert args.threads > 0
    return args


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

# Both query and target intervals are 0-based, start inclusive, and end explusive
# E.g. the interval 0-10 is 10bp long, includes the base at index 0 but not the base at index 10
def get_intervals(aln, max_del_size=20):
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

    # aln.reference_start is 0-indexed 
    tstart = aln.reference_start

    intervals = list() # list of exonic intervals of the alignment
    qstart_c = qstart # current interval's start on query
    qend_c = qstart # current interval's end on query
    tstart_c = tstart  # current interval's start on target
    tend_c = tstart # current interval's end on target
    interval_cigar = list() # current interval's list of cigar operations 
    for t, c in cigar:
        assert(0 <= t < 10), t
        # Treat any deletion (cigar D) longer than max_del_size as a target skip (cigar N)
        if t == pysam.CDEL and c > max_del_size:
            t = pysam.CREF_SKIP
        if t in exon_consuming:
            interval_cigar.append((t, c))
        if t == pysam.CDEL:
            tend_c += c
        elif t == pysam.CINS:
            qend_c += c
        elif t in target_and_query_consuming:
            tend_c += c
            qend_c += c
        # End of the current interval
        if t == pysam.CREF_SKIP:
            intervals.append((
                tstart_c,
                tend_c,
                qstart_c,
                qend_c,
                interval_cigar,
            ))
            assert(sum(c for t,c in interval_cigar if t in query_consuming) == qend_c-qstart_c)
            assert(sum(c for t,c in interval_cigar if t in target_consuming) == tend_c-tstart_c)
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
        assert(sum(c for t,c in interval_cigar if t in query_consuming) == qend_c-qstart_c)
        assert(sum(c for t,c in interval_cigar if t in target_consuming) == tend_c-tstart_c)
    # X = sum(c for t, c in cigar if t in target_consuming)
    # Y = sum(e-s for s,e,_,_,_ in intervals)
    # if X!=Y:
    #     print(X,Y,X-Y)
    #     print(aln.cigarstring)
    #     print(aln.query_name)
    #     for idx,(ts,te,qs,qe,C) in enumerate(intervals):
    #         print(idx,te-ts,sum(c for t,c in C if t in target_consuming))
    #         print(idx,qe-qs,sum(c for t,c in C if t in query_consuming))
    #     assert False
    return intervals
    # return list(fix_intervals(intervals))


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

def break_tint(tint, reads):
    rids = tint['rids']
    intervals = tint['intervals']
    start = intervals[0][0]
    end = intervals[-1][1]
    pos_to_intrv = np.zeros(end-start, dtype=int)
    pos_to_intrv[:] = len(intervals)
    intrv_to_rids = [set() for _ in intervals]
    rid_to_intrvs = {rid:set() for rid in rids}
    for idx,(s,e) in enumerate(intervals):
        pos_to_intrv[s-start:e-start] = idx
    edges = dict()
    for rid in rids:
        read = reads[rid]
        alns = read['intervals']
        for aln in alns:
            s = aln[0]
            v1 = pos_to_intrv[s-start]
            intrv_to_rids[v1].add(rid)
            rid_to_intrvs[rid].add(v1)
        for a1,a2 in zip(alns[:-1],alns[1:]):
            junc_start = a1[1]
            junc_end = a2[0]
            v1 = pos_to_intrv[junc_start-start-1]
            v2 = pos_to_intrv[junc_end-start]
            assert v1<=v2<len(intervals), (
                junc_start,
                junc_end,
                v1,
                v2,
            )
            edges[(v1,v2)] = edges.get((v1,v2), 0) + 1
    
    graph_edges = [(u,v) for (u,v),w in edges.items() if w >=2]
    G = nx.Graph()
    G.add_nodes_from(range(len(intervals)))
    G.add_edges_from(graph_edges)
    for c in nx.connected_components(G):
        c_rids = set()
        for i in c:
            c_rids.update(intrv_to_rids[i])
        if len(c_rids) > 2:
            rid_intrvs = set()
            for rid in c_rids:
                rid_intrvs.update(rid_to_intrvs[rid])
            rid_intrvs = sorted(rid_intrvs)
            yield dict(
                intervals=[intervals[i] for i in rid_intrvs],
                rids=sorted(c_rids),
            )

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
        # for rid in rids:
        #     reads[rid]['tint'] = len(multi_tints)
        multi_tints.append(dict(intervals=sorted(
            group_intervals), rids=sorted(rids)))
    assert all(enqueued)
    fin_multi_tints = list()
    for tint in multi_tints:
        if len(tint['intervals']) < 100 and len(tint['rids']) < 1500:
            fin_multi_tints.append(tint)
        else:
            X = 0
            for t in break_tint(tint,reads):
                X+=1
                fin_multi_tints.append(t)
    return fin_multi_tints


def split_reads(read_files, rname_to_tint, contigs, outdir, threads):
    outfiles = {c: open('{}/{}/reads.tsv'.format(outdir, c), 'w+')
                for c in contigs}
    for read_file in read_files:
        print('[freddie_split] Splitting reads:', read_file)
        if read_file.endswith('.gz'):
            read_file_open = gzip.open(read_file, 'rt')
        else:
            read_file_open = open(read_file, 'r')
        for idx, line in enumerate(read_file_open):
            if idx == 0:
                if line[0] == '@':
                    mod = 4
                elif line[0] == '>':
                    mod = 2
                else:
                    assert False, 'Invalid fasta/q file ' + read_file
            if idx % mod == 0:
                rname = line.rstrip().split()[0][1:]
                if rname in rname_to_tint:
                    contig = rname_to_tint[rname]['contig']
                    rid = rname_to_tint[rname]['rid']
                    tint_ids = rname_to_tint[rname]['tint_ids']
                else:
                    contig = None
            if idx % mod == 1 and contig != None:
                seq = line.rstrip()
                for tint_id in tint_ids:
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
        print('[freddie_split] Sorting contig {}...'.format(c))
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
                print('[freddie_split] Splitting contig {}'.format(contig))
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
        if not read['name'] in rname_to_tint:
            rname_to_tint[read['name']] = dict(
                contig = contig,
                rid = rid,
                tint_ids = list()
            )
        assert rname_to_tint[read['name']]['contig'] == contig
        assert rname_to_tint[read['name']]['rid'] == rid, (contig, rid, read['name'], rname_to_tint[read['name']]['rid'])
        rname_to_tint[read['name']]['tint_ids'].append(tint_id)
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

    contigs = {
        x['SN']
        for x in pysam.AlignmentFile(args.bam, 'rb').header['SQ']
        if x['LN'] > args.contig_min_size 
    }
    assert len(contigs) > 0, f'No contigs are left! Try checking BAM header or --contig-min-size parameter'
    args.threads = min(args.threads, len(contigs))
    split_args = list()
    for contig in sorted(contigs, reverse=True):
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
                rname_to_tint=rname_to_tint, contigs=final_contigs, outdir=args.outdir, threads=args.threads)

if __name__ == "__main__":
    main()
