#!/usr/bin/env python3
import argparse
import itertools
from Bio.Seq import Seq

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run NanoSim with normal or uniform distribution for the sequenced transcripts")
    parser.add_argument("-tt",
                        "--transcript-tsv",
                        type=str,
                        required=True,
                        help="Path to transcripts.tsv")
    parser.add_argument("-tf",
                        "--transcript-fasta",
                        type=str,
                        required=True,
                        help="Path to transcripts.fasta")
    parser.add_argument("-ns",
                        "--nanosim",
                        type=str,
                        required=True,
                        help="Path to NanoSim simulator.py")
    parser.add_argument("-tr",
                        "--nanosim-train-prefix",
                        type=str,
                        required=True,
                        help="Path to NanoSim train files prefix (e.g. ./train)")
    parser.add_argument("-d",
                        "--intermediate-directory",
                        type=str,
                        required=True,
                        help="Intermediate directory where NanoSim simulated files per transcript will be stored")
    parser.add_argument("-c",
                        "--read-count",
                        type=int,
                        default=20,
                        help="Number of reads to generate in total")
    parser.add_argument("-f",
                        "--read-distribution",
                        choices=['normal', 'uniform'],
                        type=str,
                        default='normal',
                        help="Number of reads to generate in total")
    parser.add_argument("-or",
                        "--output-oriented-reads",
                        type=str,
                        required=True,
                        help="Output FASTQ file for NanoSim reads of all transcripts oriented")
    parser.add_argument("-ot",
                        "--output-oriented-tsv",
                        type=str,
                        required=True,
                        help="Output file for the TSV file of the oriented reads")
    args = parser.parse_args()
    return args

def get_transcript_infos(transcript_tsv, transcript_fasta):
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
    fasta = open(transcript_fasta)
    while True:
        line = fasta.readline()
        if not len(line) > 0:
            break
        tid = line.strip().split()[0].lstrip('>')
        line = fasta.readline()
        transcript_infos[tid]['seq'] = line.strip()
    return transcript_infos

def run_nanosim(transcript_infos, nanosim, nanosim_train_prefix, intermediate_directory):
    from os import makedirs
    import subprocess
    makedirs(intermediate_directory, exist_ok=True)
    for tid in transcript_infos:
        if transcript_infos[tid]['expression_count'] < 1:
            continue
        temp_fasta_path = '{}/{}.fasta'.format(intermediate_directory, tid.replace(' ', '_'))
        temp_fasta = open(temp_fasta_path, 'w+')
        print('>{}'.format(tid), file=temp_fasta)
        print(transcript_infos[tid]['seq'], file=temp_fasta)
        temp_fasta.close()
        cmd = '{nanosim} linear -r {fasta} -c {train_prefix} -o {out_prefix} -n {read_count}'.format(
            nanosim=nanosim,
            fasta=temp_fasta_path,
            train_prefix=nanosim_train_prefix,
            out_prefix='{}/{}_simulated'.format(intermediate_directory, tid.replace(' ', '_')),
            read_count=transcript_infos[tid]['expression_count']
        )
        print(cmd)
        subprocess.run(cmd.split())

# Modified from https://www.geeksforgeeks.org/python-make-a-list-of-intervals-with-sequential-numbers/
def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    func = lambda t: t[1] - t[0]
    for key, group in itertools.groupby(enumerate(iterable),func):
        group = list(group)
        yield (group[0][1], group[-1][1]+1)

def orient_and_merge_nanosim(transcript_infos, intermediate_directory, out_oriented_reads, out_tsv):
    out_oriented_reads = open(out_oriented_reads, 'w+')
    out_tsv = open(out_tsv, 'w+')
    for tid in transcript_infos:
        if transcript_infos[tid]['expression_count'] < 1:
            continue
        fasta = open('{}/{}_simulated_reads.fasta'.format(intermediate_directory, tid.replace(' ', '_')))
        while True:
            line = fasta.readline().rstrip()
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
            aligned = line[2]
            rid = int(line[3])
            strand = line[4]
            ssc = int(line[5])
            aln = int(line[6])
            esc = int(line[7])

            oriented = 'unoriented'
            seq = fasta.readline().rstrip()
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


def generate_transcript_counts(transcript_infos, read_count, read_distribution):
    if read_distribution == 'normal':
        for tid in transcript_infos:
            transcript_infos[tid]['expression_count'] = 0
        from random import choice
        tids = list(transcript_infos.keys())
        for _ in range(read_count):
            transcript_infos[choice(tids)]['expression_count'] += 1
    elif read_distribution == 'normal':
        from math import floor
        quotient = floor(read_count/len(transcript_infos))
        remainder = read_count%len(transcript_infos)
        for tid in transcript_infos:
            transcript_infos[tid] = quotient
            transcript_infos[tid] += (remainder > 0)
            remainder -= 1
            print(tid, transcript_infos[tid].keys())

def main():
    args = parse_args()
    print(args)
    transcript_infos = get_transcript_infos(
        transcript_tsv=args.transcript_tsv,
        transcript_fasta=args.transcript_fasta
    )
    generate_transcript_counts(
        transcript_infos=transcript_infos,
        read_count=args.read_count,
        read_distribution=args.read_distribution
    )
    run_nanosim(
        transcript_infos=transcript_infos,
        nanosim=args.nanosim,
        nanosim_train_prefix=args.nanosim_train_prefix,
        intermediate_directory=args.intermediate_directory.rstrip('/')
    )
    orient_and_merge_nanosim(
        transcript_infos=transcript_infos,
        intermediate_directory=args.intermediate_directory.rstrip('/'),
        out_oriented_reads=args.output_oriented_reads,
        out_tsv=args.output_oriented_tsv
    )

if __name__ == "__main__":
    main()
