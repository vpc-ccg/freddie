#!/usr/bin/env python3
import argparse
import math
import statistics

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import upsetplot

def parse_args():
    parser = argparse.ArgumentParser(
        description="Outputs Upset plots for different tools isoform predictions")
    parser.add_argument("-a",
                        "--gtfs",
                        type=str,
                        nargs='+',
                        required=True,
                        help="Space separated list of annotaion GTF files")
    parser.add_argument("-n",
                        "--names",
                        type=str,
                        nargs='+',
                        required=True,
                        help="Space separated list of names for the annotaion GTF files")
    parser.add_argument("-t",
                        "--title",
                        type=str,
                        default='',
                        help="Title of the plot. Default: ''")
    parser.add_argument("-j",
                        "--jump",
                        type=int,
                        default=10,
                        help="The distance in bases for merging two splice sites. Default: 10")
    parser.add_argument("-x",
                        "--hides",
                        nargs='+',
                        default=list(),
                        help="Space separated list of names to hide if appearing alone. Default: None")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default='freddie_venns',
                        help="Output files prefix. Default: freddie_venns")
    parser.add_argument("-T",
                        "--truth",
                        type=str,
                        required=True,
                        help='Name associated with ground truth GTF to be used for accuracy stats.'
    )
    args = parser.parse_args()
    assert len(args.names) == len(args.gtfs)
    assert all(x in args.names for x in args.hides)
    args.hides = {(x,) for x in args.hides}
    return args


def get_juncs(gtfs):
    starts = dict()
    ends = dict()
    for tool,gtf in gtfs.items():
        for line in open(gtf):
            if line[0]=='#':
                continue
            line = line.rstrip().split('\t')
            if line[2]!='exon':
                continue
            chrom = line[0]
            if not chrom in starts:
                starts[chrom] = set()
                ends[chrom] = set()
            starts[chrom].add(int(line[3]))
            ends[chrom].add(int(line[4]))
    for chrom in starts:
        starts[chrom] = sorted(starts[chrom])
        ends[chrom] = sorted(ends[chrom])
    return starts,ends

def interval_extract(l, jump=1):
    idx = 0
    while (idx < len(l)):
        group = list()
        while True:
            group.append(l[idx])
            if idx == len(l)-1: # END of the list
                break
            if l[idx]+jump < l[idx + 1]: # next item is too far
                break
            idx += 1
        yield group
        idx += 1

def get_junc_maps(starts, ends, jump):
    starts_map = dict()
    ends_map = dict()

    for chrom in starts.keys():
        starts_map[chrom] = dict()
        ends_map[chrom] = dict()
        for group in interval_extract(starts[chrom], jump=jump):
            target = math.floor(statistics.median(group))
            for pos in group:
                starts_map[chrom][pos] = target
        for group in interval_extract(ends[chrom], jump=jump):
            target = math.floor(statistics.median(group))
            for pos in group:
                ends_map[chrom][pos] = target
    return starts_map,ends_map


def get_tool_isoforms(gtfs,starts_map,ends_map):
    tool_isoforms = dict()
    for tool,gtf in gtfs.items():
        print(f'Reading {tool}...', end=' ')
        tool_isoforms[tool] = dict()
        isoforms = tool_isoforms[tool]
        for line in open(gtf):
            if line[0]=='#':
                continue
            line = line.rstrip().split('\t')
            if line[2]!='exon':
                continue
            chrom = line[0]
            start = starts_map[chrom][int(line[3])]
            end = ends_map[chrom][int(line[4])]
            info = line[8]
            info = [x.strip().split(' ') for x in info.strip(';').split(';')]
            info = {x[0]:x[1].strip('"') for x in info}
            tid = info['transcript_id']
            
            if not tid in isoforms:
                isoforms[tid] = dict(
                    chrom = chrom,
                    exon_intervals = list(),
                )
            isoforms[tid]['exon_intervals'].append((start,end))
        for tid,isoform in isoforms.items():
            chrom = isoforms[tid]['chrom']
            exons = tuple(sorted(isoforms[tid]['exon_intervals']))
            isoforms[tid]['exons'] = (chrom, exons)
            introns = list()
            for (s1,e1),(s2,e2) in zip(exons[:-1],exons[1:]):
                introns.append((e1,s2))
            introns = tuple(introns)
            isoforms[tid]['introns'] = (chrom, introns)
        print('done')

    return tool_isoforms

def get_all_isoforms(tool_isoforms):
    all_isoforms_exonic = dict()
    all_isoforms_intronic = dict()
    for name,isoforms in tool_isoforms.items():
        for tid,isoform in isoforms.items():
            exons = isoform['exons']
            introns = isoform['introns']
            if not exons in all_isoforms_exonic:
                all_isoforms_exonic[exons] = {n:0 for n in tool_isoforms.keys()}
            all_isoforms_exonic[exons][name] += 1
            if not introns in all_isoforms_intronic:
                all_isoforms_intronic[introns] = {n:0 for n in tool_isoforms.keys()}
            all_isoforms_intronic[introns][name] += 1
    return all_isoforms_exonic,all_isoforms_intronic



def get_venns(tool_isoforms):
    feat_type_to_venn = {
        'exons':dict(),
        'introns':dict(),    
    }
    for feat_type,venn in feat_type_to_venn.items():
        isoform_key_to_tool_counts = dict()
        for tool,isoforms in tool_isoforms.items():
            for tid,isoform in isoforms.items():
                isoform_key = isoform[feat_type]
                if not isoform_key in isoform_key_to_tool_counts:
                    isoform_key_to_tool_counts[isoform_key] = {t:0 for t in tool_isoforms.keys()}
                isoform_key_to_tool_counts[isoform_key][tool] += 1
        for isoform_key,tool_counts in isoform_key_to_tool_counts.items():
            if isoform_key == tuple():
                continue
            venn_key = tuple(t for t,c in tool_counts.items() if c > 0)
            venn[venn_key] = venn.get(venn_key, 0) + 1
    return feat_type_to_venn

def upset_plot_venns(feat_type_to_venn, prefix, hides):
    for feat_type,venn in feat_type_to_venn.items():
        venn = {k:v for k,v in venn.items() if not k in hides}
        fig = plt.figure(figsize=(10,30), )
        df = upsetplot.from_memberships(
            [list(x[0]) for x in sorted(venn.items())],
            [x[1] for x in sorted(venn.items())],
        )
        upsetplot.plot(df, show_counts='%d', fig=fig)
        plt.savefig('{}.{}.pdf'.format(prefix, feat_type))

def output_accuracy_stats(gtfs, feat_type_to_venn, truth_name, prefix):
    for feat_type,venn in feat_type_to_venn.items():
        outfile = open(f'{prefix}.{feat_type}.tsv', 'w+')
        record = [
            'Tool', 
            'F1-score', 
            'Precision', 
            'Recall', 
            'True isoforms', 
            'False isoforms',
        ]
        outfile.write('\t'.join(record))
        outfile.write('\n')
        truth_count = sum(c for x,c in venn.items() if truth_name in x)
        for t in gtfs:
            d = {False:0, True:0}
            for x,c in venn.items():
                if t in x:
                    d[truth_name in x] += c
            T = d[True]
            F = d[False]
            p = T/(T+F)
            r = T/(truth_count)
            f1 = 2*p*r/(p+r)
            record = [
                f'{t:}',
                f'{f1:2.2%}',
                f'{p:2.2%}',
                f'{r:2.2%}',
                f'{T}',
                f'{F}',
            ]
            outfile.write('\t'.join(record))
            outfile.write('\n')
        outfile.close()

def main():
    args = parse_args()
    args.gtfs = {
        t:g for t,g in zip(args.names, args.gtfs)
    }
    starts,ends = get_juncs(args.gtfs)
    print(f'Starts: {sum(len(v) for v in starts.values())}')
    print(f'Ends: {sum(len(v) for v in ends.values())}')
    starts_map, ends_map = get_junc_maps(starts, ends, args.jump)    
    tool_isoforms = get_tool_isoforms(args.gtfs, starts_map, ends_map)
    feat_type_to_venn = get_venns(tool_isoforms)
    upset_plot_venns(feat_type_to_venn, args.output, args.hides)
    output_accuracy_stats(args.gtfs, feat_type_to_venn, args.truth, args.output)


if __name__ == "__main__":
    main()
