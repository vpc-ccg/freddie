#!/usr/bin/env python3
import argparse
import math
import statistics
import os
from collections import Counter

import upsetplot
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter,PercentFormatter
import pandas as pd

mpl.use('Agg')


def parse_args():
    parser = argparse.ArgumentParser(
        description="Outputs Upset plots for different tools isoform predictions with short-read validation")
    parser.add_argument("-a",
                        "--gtfs",
                        type=str,
                        nargs='+',
                        required=True,
                        help="Space separated list of annotaion GTF files")
    parser.add_argument("-sj",
                        "--short-juncs",
                        type=str,
                        nargs='+',
                        required=True,
                        help="Space separated list of STAR mapper SJ.out.tab files")
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
    args = parser.parse_args()
    assert len(args.names) == len(args.gtfs)
    args.gtfs = dict(zip(
        args.names,
        args.gtfs
    ))
    assert all(x in args.names for x in args.hides)
    args.hides = {(x,) for x in args.hides}
    return args


def get_validated_junc(junc_files, jump):
    juncs_valid = dict()
    for junc_file in junc_files:
        for l in open(junc_file):
            chrom,s,e = l.rstrip().split('\t')[:3]
            if not chrom in juncs_valid:
                juncs_valid[chrom] = set()
            s = int(s)
            e = int(e)
            for i in range(s-jump,s+jump):
                for j in range(e-jump,e+jump):
                    juncs_valid[chrom].add((i,j))
    return juncs_valid

def get_juncs(tool_isoforms):
    starts = dict()
    ends = dict()
    for tool,isoforms in tool_isoforms.items():
        for isoform in isoforms.values():
            chrom = isoform['chrom']
            if not chrom in starts:
                starts[chrom] = set()
                ends[chrom] = set()
            for s,e in isoform['exons']:
                starts[chrom].add(s)
                ends[chrom].add(e)
    for chrom in starts:
        starts[chrom] = sorted(starts[chrom])
        ends[chrom] = sorted(ends[chrom])
    print(f'Starts: {sum(len(v) for v in starts.values())}')
    print(f'Ends: {sum(len(v) for v in ends.values())}')
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


def get_tools_isoforms(gtfs):
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
            start = int(line[3])
            end = int(line[4])
            info = line[8]
            info = [x.strip().split(' ') for x in info.strip(';').split(';')]
            info = {x[0]:x[1].strip('"') for x in info}
            tid = info['transcript_id']
            if not tid in isoforms:
                isoforms[tid] = dict(
                    chrom = chrom,
                    exons = list(),
                )
            isoforms[tid]['exons'].append((start,end))
        print('done')

    return tool_isoforms

def get_canonized_isoforms(tool_isoforms, starts_map, ends_map):
    canonized_isoforms_e = dict()
    canonized_isoforms_i = dict()
    for tool,isoforms in tool_isoforms.items():
        print('Canonizing', tool)
        for tid,isoform in isoforms.items():
            chrom = isoform['chrom']
            canonized_exons = tuple(sorted([(starts_map[chrom][s], ends_map[chrom][e]) for s,e in isoform['exons']]))
            if not (chrom,canonized_exons) in canonized_isoforms_e:
                canonized_isoforms_e[(chrom,canonized_exons)] = dict(
                    chrom = chrom,
                    exons = canonized_exons,
                    tool_isoforms = set()
                )
            canonized_isoforms_e[(chrom,canonized_exons)]['tool_isoforms'].add((tool,tid))
            
            canonized_introns = tuple(
                [(e1,s2) for (_,e1),(s2,_) in zip(canonized_exons[:-1],canonized_exons[1:])]
            )
            if len(canonized_introns) == 0:
                continue
            if not (chrom,canonized_introns) in canonized_isoforms_i:
                canonized_isoforms_i[(chrom,canonized_introns)] = dict(
                    chrom = chrom,
                    introns = canonized_introns,
                    tool_isoforms = set()
                )
            canonized_isoforms_i[(chrom,canonized_introns)]['tool_isoforms'].add((tool,tid))
    return canonized_isoforms_e,canonized_isoforms_i

def get_venns(tool_isoforms):
    feat_type_to_venn = {
        'exons':dict(),
        'introns':dict(),    
    }
    for feat_type,venn in feat_type_to_venn.items():
        isoform_key_to_tool_counts = dict()
        for tool,isoforms in tool_isoforms.items():
            for tid,isoform in isoforms.items():
                isoform_key = tuple(isoform[feat_type])
                if not isoform_key in isoform_key_to_tool_counts:
                    isoform_key_to_tool_counts[isoform_key] = {t:0 for t in tool_isoforms.keys()}
                isoform_key_to_tool_counts[isoform_key][tool] += 1
        for isoform_key,tool_counts in isoform_key_to_tool_counts.items():
            if isoform_key == tuple():
                continue
            venn_key = tuple(t for t,c in tool_counts.items() if c > 0)
            venn[venn_key] = venn.get(venn_key, 0) + 1
    return feat_type_to_venn

def get_venn_from_canonized(isoforms, validation_k='introns', validations=dict()):
    venn = dict()
    venn_v = dict()
    for i in isoforms.values():
        try:
            v = all(se in validations[i['chrom']] for se in i[validation_k])
        except:
            v = False
        k = tuple(sorted({t for t,_ in i['tool_isoforms']}))
        if not k in venn:
            venn[k] = 0
            venn_v[k] = 0
        venn[k] += 1
        if v:
            venn_v[k]+=1
    return venn,venn_v

def get_canonized_introns(tool_isoforms, starts_map, ends_map):
    canonized_introns = dict()
    for tool,isoforms in tool_isoforms.items():
        print('Canonizing', tool)
        for tid,isoform in isoforms.items():
            chrom = isoform['chrom']
            canonized_exons = tuple(sorted([(starts_map[chrom][s], ends_map[chrom][e]) for s,e in isoform['exons']]))
            for (_,e1),(s2,_) in zip(canonized_exons[:-1],canonized_exons[1:]):
                k = (chrom,(e1,s2)) 
                if not k in canonized_introns:
                    canonized_introns[k] = dict(
                        chrom=chrom,
                        introns=[(e1,s2)],
                        tool_isoforms = set(),
                    )
                canonized_introns[k]['tool_isoforms'].add((tool,tid))
    return canonized_introns

def upset_plot_venn(venn, venn_v, prefix, hides=set(), title=''):
    venn = {k:v for k,v in venn.items() if not k in hides}
    venn_v_d = {venn[k]:venn_v[k] for k in venn}
    
    c_total = Counter()
    c_valid = Counter()
    for k in venn.keys():
        all_c = venn[k]
        val_c = venn_v[k]
        for t in k:
            c_total[t]+=all_c
            c_valid[t]+=val_c
    venn_totals_v_d = {c_total[t]:c_valid[t] for t in c_total.keys()} 

    fig = plt.figure(figsize=(10,35), )
    fig.suptitle(title)
    df = upsetplot.from_memberships(
        [list(x[0]) for x in sorted(venn.items())],
        [x[1] for x in sorted(venn.items())],
    )
    x = upsetplot.plot(df, show_counts='%d', fig=fig)
    for c in x['intersections'].get_children()[len(venn):len(venn)*2]:
        c.set_rotation(30)
    for c in x['intersections'].get_children()[:len(venn)]:
        h = venn_v_d.get(c.get_height(), 0)
        if h == 0:
            continue
        x['intersections'].add_patch(mpl.patches.Rectangle(
            xy=c.get_xy(),
            width=c.get_width(),
            height=h,
            color='green',
            zorder=100000,
        ))
        x['intersections'].text(
            x = c.get_xy()[0]-.5,
            y = -max(venn.values())*.135,
            s = h,
            c = 'green',
            zorder=100001,
            rotation=30,
        )    
    for c in x['totals'].get_children()[:len(venn_totals_v_d)]:
        w = venn_totals_v_d.get(c.get_width(),0)
        if w == 0:
            continue
        x['totals'].add_patch(mpl.patches.Rectangle(
            xy=c.get_xy(),
            width=w,
            height=c.get_height(),
            color='green',
            zorder=100000,
        ))
        x['totals'].text(
            x = c.get_xy()[0]+0,
            y = c.get_xy()[1]+.55,
            s = w,
            c = 'green',
            zorder=100001,
            rotation=0,
            horizontalalignment='right',
        )
    plt.savefig('{}.pdf'.format(prefix))

def output_validation_rate(gtfs, venn_i, venn_i_v, venn_introns, venn_introns_v, outpath):
    outfile = open(outpath, 'w+')
    total_to_valid = dict()
    for venn,venn_v,title in [
        (venn_i,venn_i_v,'Isoforms'),
        (venn_introns,venn_introns_v,'Introns'),
    ]:
        print(title, file=outfile)
        record = [
            'Tool'.ljust(10),
            'Total',
            'Valid',
            '%',
        ]
        print('\t'.join(record), file=outfile)
        for t in gtfs:
            total = 0
            validated = 0
            for v in venn:
                if not any(t in x for x in v):
                    continue
                total+=venn[v]
                validated+=venn_v[v]
            record = [
                t,
                str(total),
                str(validated),
                f'{validated/total:2.2%}',
            ]
            total_to_valid[total] = validated

            print('\t'.join(record), file=outfile)
    outfile.close()

def main():
    args = parse_args()
    tool_isoforms = get_tools_isoforms(args.gtfs)
    juncs_valid = get_validated_junc(junc_files=args.short_juncs, jump=args.jump)
    starts,ends = get_juncs(tool_isoforms)
    starts_map, ends_map = get_junc_maps(starts, ends, args.jump)

    canonized_isoforms_e,canonized_isoforms_i = get_canonized_isoforms(
        tool_isoforms,
        starts_map,
        ends_map
    )
    canonized_introns = get_canonized_introns(
        tool_isoforms,
        starts_map,
        ends_map
    )
    
    venn_i,venn_i_v = get_venn_from_canonized(canonized_isoforms_i, 'introns', juncs_valid)
    print('Plotting isoform_intronic')
    upset_plot_venn(
        venn     = venn_i,
        venn_v   = venn_i_v,
        prefix   = f'{args.output}.isoform_intronic',
        hides    = args.hides,
        title    = f'Real data isoforms validation with short reads', 
    )

    print('Plotting introns')
    venn_introns,venn_introns_v = get_venn_from_canonized(canonized_introns, 'introns', juncs_valid)
    upset_plot_venn(
        venn   = venn_introns,
        venn_v = venn_introns_v,
        prefix = f'{args.output}.introns',
        hides  = args.hides,
        title  = 'Real data intron validation with short reads'
    )

    print('Outputting validation rates')
    output_validation_rate(
        gtfs = args.gtfs,
        venn_i = venn_i,
        venn_i_v = venn_i_v,
        venn_introns = venn_introns,
        venn_introns_v = venn_introns_v,
        outpath = f'{args.output}.validation.tsv'
    )


if __name__ == "__main__":
    main()
