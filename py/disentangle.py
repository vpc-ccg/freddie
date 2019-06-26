#!/usr/bin/env python3
import argparse
import random
from statistics import mean
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    parser.add_argument("-t",
                        "--transcript-tsv",
                        type=str,
                        required=True,
                        help="Path to TSV file of annotation transcripts")
    parser.add_argument("-e",
                        "--total-expression",
                        type=int,
                        default=500,
                        help="How many full length transcript reads to simulate in total")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        help="Output TXT file")
    args = parser.parse_args()
    return args

def plot_coverage(coverage, ticks=list(), outpath='out.pdf'):
    plt.figure(figsize=(30,4))
    # plt.plot(range(len(coverage)), coverage)
    plt.scatter(range(len(coverage)), coverage)
    height = max(coverage)
    for tick in ticks:
        plt.plot([tick,tick], [0,height], 'k--', alpha=0.2)
    plt.tight_layout()
    plt.savefig(outpath)

def main():
    random.seed(42)
    args = parse_args()
    out_file = open(args.output, 'w+')

    sites = set()
    transcripts = list()

    for line in open(args.transcript_tsv):
        line = line.rstrip().split('\t')
        exons = [(int(e.split('-')[0]),int(e.split('-')[1])) for e in line[3].split(',')]
        for exon in exons:
            sites.add(exon[0])
            sites.add(exon[1])
        transcripts.append(exons)
    sites.add(0)
    sites = sorted(sites)
    g_len = sites[-1] + 10
    sites.append(g_len)

    tids = [i for i in range(len(transcripts))]
    coverage = [0 for _ in range(g_len)]

    picked_transcripts = random.choices(tids, k=args.total_expression)
    for tid in picked_transcripts:
        for exon in transcripts[tid]:
            for i in range(exon[0], exon[1]+1):
                coverage[i] += 1
    plot_coverage(coverage=coverage, ticks=sites, outpath=args.output+'.pdf')
    canonical_exons = list()
    for i in range(1, len(sites)):
        start = sites[i-1]+1
        end = sites[i]
        if start == end or end > len(coverage):
            continue
        canonical_exons.append(mean(coverage[start:end]))
        print(start, end, canonical_exons[-1])
    plot_coverage(coverage=canonical_exons, outpath=args.output+'.2.pdf')



    out_file.close()


if __name__ == "__main__":
    main()
