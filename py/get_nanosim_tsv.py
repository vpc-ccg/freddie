#!/usr/bin/env python3
import argparse
import itertools
from Bio.Seq import Seq

def parse_args():
    parser = argparse.ArgumentParser(
        description="Get oriented version of nanosim reads their TSV file ")
    parser.add_argument("-nsr",
                        "--nanosim-reads",
                        type=str,
                        required=True,
                        help="Path to NanoSim reads in FASTQ format")
    parser.add_argument("-t",
                        "--transcript-tsv",
                        type=str,
                        required=True,
                        help="Transcripts TSV file")
    parser.add_argument("-or",
                        "--oriented-reads",
                        type=str,
                        required=True,
                        help="Output file for oriented reads")
    parser.add_argument("-ot",
                        "--oriented-tsv",
                        type=str,
                        required=True,
                        help="Output file for the TSV file")
    args = parser.parse_args()
    return args

def get_transcript_infos(transcript_tsv):
    transcript_infos = dict()
    for line in open(transcript_tsv):
        line = line.rstrip().split('\t')
        tid = line[0]
        transcript_infos[tid] = dict()
        transcript_infos[tid]['chr']    = line[1]
        transcript_infos[tid]['strand'] = line[2]
        transcript_infos[tid]['exons'] = [(int(x.split('-')[0]),int(x.split('-')[1])) for x in line[3].rstrip(',').split(',')]
        transcript_infos[tid]['length'] = sum([x[1]-x[0] for x in transcript_infos[tid]['exons']])
        transcript_infos[tid]['genic_pos'] = [0]*transcript_infos[tid]['length']
        transcript_pos = 0
        for (start,end) in transcript_infos[tid]['exons']:
            for genic_pos in range(start,end):
                transcript_infos[tid]['genic_pos'][transcript_pos] = genic_pos
                transcript_pos+=1
    return transcript_infos

# Modified from https://www.geeksforgeeks.org/python-make-a-list-of-intervals-with-sequential-numbers/
def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    func = lambda t: t[1] - t[0]
    for key, group in itertools.groupby(enumerate(iterable),func):
        group = list(group)
        yield (group[0][1], group[-1][1]+1)

def output_nanosim_reads_tsv(transcript_infos, nanosim_reads_fasta, out_oriented_reads, out_tsv):
    out_oriented_reads = open(out_oriented_reads, 'w+')
    out_tsv = open(out_tsv, 'w+')
    fasta = open(nanosim_reads_fasta)
    while True:
        line = fasta.readline()
        if len(line) == 0:
            break
        rname = line[1:].rstrip().split()[0]
        try:
            comment = line[1:].rstrip().split()[1]
        except:
            comment = ''
        line = rname.split('_')
        tid = line[0]
        start = int(line[1])
        aligned = line[2] == 'aligned'
        rid = int(line[3])
        strand = line[4]
        ssc = int(line[5])
        aln = int(line[6])
        esc = int(line[7])

        oriented = 'unoriented'
        seq = fasta.readline()
        if strand == 'R':
            oriented = 'oriented'
            ssc,esc = esc,ssc
            seq=str(Seq(seq).reverse_complement())
        rname = '{tid}_{start}_{aligned}_{rid}_{strand}_{ssc}_{aln}_{esc}_{oriented}'.format(
            tid=tid,
            start=start,
            aligned=aligned,
            rid=rid,
            strand=strand,
            ssc=ssc,
            aln=aln,
            esc=esc,
            oriented=oriented,
        )
        genic_positions = (transcript_infos[tid]['genic_pos'][start:start+aln])
        exons = list(intervals_extract(genic_positions))
        exons = ','.join(('{}-{}'.format(e[0],e[1]) for e in exons))

        print('>{}\n{}'.format(rname,seq), file=out_oriented_reads)
        print('\t'.join((rname, tid, strand, exons)), file=out_tsv)
    out_oriented_reads.close()
    out_tsv.close()

def main():
    args = parse_args()
    transcript_infos = get_transcript_infos(args.transcript_tsv)
    output_nanosim_reads_tsv(transcript_infos, args.nanosim_reads, args.oriented_reads, args.oriented_tsv)

if __name__ == "__main__":
    main()
