#!/usr/bin/env python3
import argparse
import re
import pysam
from itertools import zip_longest
from collections import namedtuple

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
                        "--output",
                        type=str,
                        default='freddie_split.tsv',
                        help="Path to output file. Default: freddie_split.tsv")
    args = parser.parse_args()
    assert args.sam_format in ['sam','bam',None]
    return args

def fix_cigar(cigar):
    fixed_cigar = list()
    fixed_cigar.append(cigar[0])
    for c,t in cigar[1:]:
        last_c,last_t = fixed_cigar[-1]
        if t==last_t:
            fixed_cigar[-1]=(last_c+c,t)
            continue
        if t in ['D','N'] and last_t in ['D','N']:
            fixed_cigar[-1]=(c+last_c,'N')
            continue
        fixed_cigar.append((c,t))
    return fixed_cigar

def get_intervals(aln):
    if aln.cigarstring == None:
        return list()
    cigar = [x for x in cigar_re.findall(aln.cigarstring)]
    assert sum(len(x[0])+len(x[1]) for x in cigar)==len(aln.cigarstring), 'Errorenous cigar for aln: {}'.format(aln)
    cigar = [(int(x[0]),x[1]) for x in cigar_re.findall(aln.cigarstring)]
    qstart = 0
    if cigar[0][1] in ['S']:
        qstart += cigar[0][0]
    qlen = 0
    for c,t in cigar:
        if t in ['I','S','M','=','X']:
            qlen+=c
    assert qlen==len(aln.query_sequence)
    qend = qlen
    if cigar[-1][1] in ['S']:
        qend -= cigar[-1][0]
    assert qend>qstart
    tstart = aln.reference_start
    tend = tstart + 1

    cigar = fix_cigar(cigar)
    for c,t in cigar:
        if t in ['D','M','=','X']:
            tend+=c
    qstart_c = qstart
    qend_c   = qstart
    tstart_c = tstart
    tend_c   = tstart

    intervals = list()
    interval_cigar = list()
    for c,t in cigar:
        if t in ['D','I','X', '=', 'M']:
            interval_cigar.append('{}{}'.format(c,t))
        if t in ['D']:
            tend_c+=c
        if t in ['I']:
            qend_c+=c
        if t in ['X', '=', 'M']:
            tend_c+=c
            qend_c+=c
        if t == 'N':
            intervals.append((
                tstart_c,
                tend_c,
                qstart_c,
                qend_c,
                ''.join(interval_cigar)
            ))
            interval_cigar = list()
            tend_c+=c
            tstart_c=tend_c
            qstart_c=qend_c
    if tstart_c < tend_c:
        intervals.append((
            tstart_c,
            tend_c,
            qstart_c,
            qend_c,
            ''.join(interval_cigar)
        ))
    return intervals

def read_sam(sam, format):
    if format=='sam':
        sam = pysam.AlignmentFile(sam, 'r')
    elif format=='bam':
        sam = pysam.AlignmentFile(sam, 'rb')

    contigs = {x['SN']:int(x['LN']) for x in sam.header['SQ']}
    reads = list()
    for aln in sam.fetch():
        if aln.is_secondary or aln.reference_name == None:
            continue
        assert aln.reference_name in contigs, '{} : {}'.format(aln.reference_name, contigs)
        reads.append(dict(
            id        = len(reads),
            name      = aln.query_name,
            contig    = aln.reference_name,
            strand    = '-' if aln.is_reverse else '+',
            intervals = get_intervals(aln),
            tint      = float('inf'),
        ))
    return contigs,reads

def get_transcriptional_intervals(reads, contig):
    intervals = list()
    start = None
    end = None
    rids = list()
    for s,e,rid in sorted((i[0],i[1],read['id']) for read in reads if read['contig']==contig for i in read['intervals']):
        if (start,end) == (None,None):
            start,end = s,e
        if s >= end:
            intervals.append(dict(
                tint=-1,
                start=start,
                end=end,
                rids=rids,
            ))
            start = s
            end = e
            rids = list()
        assert start<=s
        end = max(end, e)
        rids.append(rid)
        reads[rid]['tint'] = min(reads[rid]['tint'],len(intervals))
    intervals.append(dict(
        tint=-1,
        start=start,
        end=end,
        rids=rids,
    ))
    final_tints = dict()
    for idx,interval in enumerate(intervals):
        tint = min(reads[rid]['tint'] for rid in interval['rids'])
        interval['tint'] = tint
        for rid in rids:
            reads[rid]['tint']=tint
        if not tint in final_tints:
            final_tints[tint]=list()
        final_tints[tint].append(idx)
    final_tints = sorted(final_tints.values())
    final_tint_to_intervals = list()
    for idx,tints in enumerate(final_tints):
        rids = sorted({rid for tint in tints for rid in intervals[tint]['rids']})
        final_intervals = [(intervals[tind]['start'],intervals[tind]['end']) for tind in tints]
        for rid in rids:
            reads[rid]['tint']=idx
        final_tint_to_intervals.append(dict(intervals=final_intervals, rids=rids))
    return final_tint_to_intervals

def main():
    args = parse_args()

    if args.sam_format == None:
        args.sam_format = 'sam'
        if args.sam.endswith('.bam'):
            args.sam_format = 'bam'
    contigs,reads=read_sam(sam=args.sam, format=args.sam_format)
    outfile = open(args.output, 'w')
    for contig in contigs.keys():
        tints = get_transcriptional_intervals(reads=reads, contig=contig)
        for idx,tint in enumerate(tints):
            record = list()
            record.append('#{}:{}'.format(contig,idx))
            record.append(','.join('{}:{}'.format(s,e) for s,e in tint['intervals']))
            record.append(str(len(tint['rids'])))
            outfile.write('\t'.join(record))
            outfile.write('\n')
    for read in sorted(reads, key=lambda x:(x['contig'],x['tint'])):
        record = list()
        record.append(str(read['id']))
        record.append(read['name'])
        record.append(read['contig'])
        record.append(str(read['tint']))
        for interval in read['intervals']:
            record.append('{}-{}:{}-{}:{}'.format(*interval))
        outfile.write('\t'.join(record))
        outfile.write('\n')
    outfile.close()

if __name__ == "__main__":
    main()
