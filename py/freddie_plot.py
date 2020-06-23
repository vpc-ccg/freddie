#!/usr/bin/env python3
from multiprocessing import Pool
import os
import argparse
import re
from operator import itemgetter

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from PyPDF2 import PdfFileMerger, PdfFileReader

colors = [
    '#a6cee3',
    '#1f78b4',
    '#b2df8a',
    '#33a02c',
    '#fb9a99',
    '#e31a1c',
    '#fdbf6f',
    '#ff7f00',
    '#cab2d6',
    '#6a3d9a',
    '#ffff99',
    '#b15928',
]
grid_width_ratios = [
    (5000, 9),
    (2000, 8),
    (1000, 7),
    (500, 6),
    (200, 5),
    (100, 4),
    (50, 3),
    (20, 2),
    (0, 1)
]

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    parser.add_argument("-a",
                        "--annotation-gtf",
                        type=str,
                        required=True,
                        help="Path to GTF file of annotations")
    parser.add_argument("-s",
                        "--segment-tsv",
                        type=str,
                        required=True,
                        help="Path to TSV file of Freddie segment")
    parser.add_argument("-c",
                        "--cluster-tsv",
                        type=str,
                        required=True,
                        help="Path to TSV file of Freddie cluster")
    parser.add_argument("-t",
                        "--threads",
                        type=int,
                        default=1,
                        help="Number of threads. Default: 1")
    parser.add_argument("-od",
                        "--out-dir",
                        type=str,
                        default='freddie_plot',
                        help="Output directory. Will be created if does not exist. Default: freddie_plot")
    args = parser.parse_args()
    args.out_dir = args.out_dir.rstrip('/')
    assert 1 <= args.threads <= 256
    return args

def plot_isoform(isoform, transcripts, plot_settings, outpath):
    grid_lens = plot_settings['grid_lens']
    segs      = plot_settings['segs']
    s_tail_order = ['SC','_poly_size','_poly_gap',]
    e_tail_order = list(reversed(s_tail_order))
    s_tail_gs_ratios = list()
    e_tail_gs_ratios = list()
    for x in s_tail_order:
        sm = 2000
        for threshold,value in grid_width_ratios:
            if sm > threshold:
                s_tail_gs_ratios.append(value)
                break
    for x in e_tail_order:
        em = 2000
        for threshold,value in grid_width_ratios:
            if em > threshold:
                e_tail_gs_ratios.append(value)
                break
    fig = plt.figure(figsize=(min(60,len(s_tail_gs_ratios)+len(grid_lens)+len(e_tail_gs_ratios)), 30), constrained_layout=False)
    out_gs = fig.add_gridspec(
        ncols=3,
        nrows=1,
        width_ratios=[sum(s_tail_gs_ratios), sum(grid_lens), sum(e_tail_gs_ratios)],
        wspace=0.05
    )
    gs = out_gs[0].subgridspec(
        ncols=len(s_tail_gs_ratios),
        nrows=2,
        height_ratios=[1,5],
        width_ratios=s_tail_gs_ratios,
        hspace=.1,
        wspace=.10,
    )
    s_axes = [fig.add_subplot(gs[1,i]) for i in range(len(s_tail_gs_ratios))]
    gs = out_gs[2].subgridspec(
        ncols=len(e_tail_gs_ratios),
        nrows=2,
        height_ratios=[1,5],
        width_ratios=e_tail_gs_ratios,
        hspace=.1,
        wspace=.10,
    )
    e_axes = [fig.add_subplot(gs[1,i]) for i in range(len(e_tail_gs_ratios))]
    gs = out_gs[1].subgridspec(
        ncols=len(grid_lens),
        nrows=2,
        width_ratios=grid_lens,
        height_ratios=[1,5],
        hspace=.1,
        wspace=0
    )
    t_axes = [fig.add_subplot(gs[0,i]) for i in range(len(grid_lens))]
    r_axes = [fig.add_subplot(gs[1,i]) for i in range(len(grid_lens))]

    ylim=len(isoform['reads'])
    for ax in [s_axes[0], s_axes[2],]:
        ax.set_xticks([50,100,500,1000,2000])
        ax.set_ylim(0,ylim)
        ax.set_yticks([])
        ax.set_xscale('log')
        ax.set_xlim(sm, 1)
    s_axes[0].set_title('Extra SC', loc='left', y=-0.025)
    s_axes[2].set_title('Gap', loc='right', y=-0.025)
    for ax in [s_axes[1],]:
        ax.set_title('Start polyA/T')
        ax.xaxis.tick_top()
        # ax.tick_params(axis='x', rotation=90)
        # ax.set_xticks([0, 10, 30, 50, 100,150])
        ax.set_ylim(0,ylim)
        ax.set_yticks([])
        ax.set_xlim(150, 0)
    for ax in [e_axes[0], e_axes[2],]:
        ax.set_xticks([50,100,500,1000,2000])
        ax.set_ylim(0,ylim)
        ax.set_yticks([])
        ax.set_xscale('log')
        ax.set_xlim(1,em)
    e_axes[0].set_title('Gap', loc='left', y=-0.025)
    e_axes[2].set_title('Extra SC', loc='right', y=-0.025)
    for ax in [e_axes[1],]:
        ax.set_title('End polyA/T')
        ax.xaxis.tick_top()
        # ax.tick_params(axis='x', rotation=90)
        # ax.set_xticks([0, 10, 30, 50, 100,150])
        ax.set_ylim(0,ylim)
        ax.set_yticks([])
        ax.set_xlim(0, 150)
    for axes,ylim in [(t_axes,len(plot_settings['plot_tids'])),(r_axes,len(isoform['reads']))]:
        for ax,(s,e) in zip(axes,segs):
            ax.grid(zorder=0)
            ax.set_ylim(0,ylim)
            ax.set_xlim(0, 1)
            ax.set_xticks([0])
            ax.set_xticklabels([s], rotation=45)
            ax.set_yticks(range(0,ylim,max(ylim//20,1)))
            if ax.is_first_col():
                pass
            elif ax.is_last_col():
                ax.yaxis.tick_right()
                ax.set_xticks([0,1])
                ax.set_xticklabels([s,e], rotation=45)
            else:
                ax.set_yticklabels([])
                # ax.set_yticks([])
    t_axes[0].set_yticklabels( [transcripts[tid]['name'] for tid in plot_settings['plot_tids']])
    t_axes[-1].set_yticklabels([transcripts[tid]['name'] for tid in plot_settings['plot_tids']])

    for p,read in enumerate(sorted(isoform['reads'], key=itemgetter('tid', 'data'))):
        if read['tid'] in plot_settings['tid_colors']:
            color = plot_settings['tid_colors'][read['tid']]
        else:
            color = 'gray'
        # for s_ax,x in zip(s_axes,s_tail_order):
        #     k = 'S{}'.format(x)
        #     if x =='_poly_size' and (strand=='+' or read['tail']['S_poly_type']=='A'):
        #         s_ax.add_patch(patches.Rectangle(
        #             xy     = (0,p),
        #             width  = read['tail'][k],
        #             height = 1,
        #             color  = 'red',
        #             hatch  = '/',
        #         ))
        #     else:
        #         s_ax.add_patch(patches.Rectangle(
        #             xy     = (0,p),
        #             width  = read['tail'][k],
        #             height = 1,
        #             color  = color,
        #         ))
        # for e_ax,x in zip(e_axes,e_tail_order):
        #     k = 'E{}'.format(x)
        #     if x =='_poly_size' and (strand=='-' or read['tail']['E_poly_type']=='T'):
        #         e_ax.add_patch(patches.Rectangle(
        #             xy     = (0,p),
        #             width  = read['tail'][k],
        #             height = 1,
        #             color  = 'red',
        #             hatch  = '/',
        #         ))
        #     else:
        #         e_ax.add_patch(patches.Rectangle(
        #             xy     = (0,p),
        #             width  = read['tail'][k],
        #             height = 1,
        #             color  = color,
        #         ))
        for aid,ax in enumerate(r_axes):
            if read['data'][aid]=='0':
                continue
            ax.add_patch(patches.Rectangle(
                xy     = (0,p),
                width  = 1,
                height = 1,
                color  = color,
            ))
    # for ax in r_axes+[sce_ax,scs_ax]:
    #     ax.hlines(y=ch_changes, xmin=0, xmax=1, alpha=0.4)
    for p,tid in enumerate(plot_settings['plot_tids']):
        for s,e in transcripts[tid]['intervals']:
            t_aid = -1
            x_0 = -1
            x_1 = -1
            for idx,(seg_s,seg_e) in enumerate(segs):
                if not (s <= seg_e and e >= seg_s):
                    continue
                t_aid = idx
                seg_l = seg_e-seg_s
                if seg_s < s:
                    x_0 = (s-seg_s)/seg_l
                else:
                    x_0 = 0.0
                if seg_e < e:
                    x_1 = 1.0
                else:
                    x_1 = (e-seg_s)/seg_l
                t_axes[t_aid].add_patch(patches.Rectangle(xy=(x_0,p),width=x_1,height=1,color=plot_settings['tid_colors'][tid]))
            # t_axes[0].text(x=-0.5, y=tid_to_p[tid]+0.5, s=tid)
    plt.savefig(outpath,bbox_inches='tight')
    print('Done with',outpath)

def get_transcripts(gtf):
    transcripts = dict()
    for line in open(gtf):
        if line[0] == '#':
            continue
        line=line.rstrip().split('\t')
        if not line[2]=='exon':
            continue
        tid = re.search(r'transcript_id \"(?P<tid>ENST\d{11})\"', line[8]).group('tid')
        interval = (int(line[3]), int(line[4]))
        if not tid in transcripts:
            transcripts[tid]=dict(
                chrom=line[0],
                intervals=list(),
                name=re.search(r'transcript_name \"(?P<tname>[\w\.-]+)\"', line[8]).group('tname')
            )
        transcripts[tid]['intervals'].append(interval)
    return transcripts

def get_tints(cluster_tsv, segment_tsv):
    tints = dict()
    rid_to_data=dict()
    for line in open(segment_tsv):
        if line[0]=='#':
            continue
        line=line.rstrip().split('\t')
        rid_to_data[int(line[0])]=line[5]
    for line in open(cluster_tsv):
        if line.startswith('#'):
            chrom,tint_id,segs=line.rstrip()[1:].split('\t')
            tint_id = int(tint_id)
            segs = segs.split(',')
            segs = [(int(s),int(e)) for s,e in zip(segs[:-1], segs[1:])]
            tints[tint_id]=dict(
                id=tint_id,
                chrom=chrom,
                segs=segs,
                isoforms={'garbage': dict(id='garbage', data='', reads=list())}
            )
        elif line.startswith('isoform_'):
            iid,tint_id,data = line.rstrip()[8:].split('\t')
            tint_id=int(tint_id)
            tints[tint_id]['isoforms'][iid]=dict(
                id=iid,
                data=data,
                reads=list(),
            )
        else:
            line=line.rstrip().split('\t')
            rid=int(line[0])
            tint=int(line[4])
            iid=line[5]
            if iid =='*':
                iid = 'garbage'
            data=rid_to_data[rid]
            tints[tint]['isoforms'][iid]['reads'].append(dict(
                rid=rid,
                name=line[1],
                tid=line[1].split('_')[0],
                chrom=line[2],
                strand=line[3],
                tint=tint,
                iid=iid,
                data=data,
                poly_tail=line[6+1+len(data):],
            ))
    for tint_id in tints.keys():
        tints[tint_id]['tids'] = set()
        for isoform in tints[tint_id]['isoforms'].values():
            for read in isoform['reads']:
                tints[tint_id]['tids'].add(read['tid'])
    return tints

def plot_tint(plot_args):
    tint,transcripts,out_dir = plot_args
    print('Outputting to',out_dir)
    os.makedirs(out_dir, exist_ok=True)
    plot_settings = dict()

    grid_lens = list()
    for s,e in tint['segs']:
        for threshold,value in grid_width_ratios:
            if e-s > threshold:
                grid_lens.append(value)
                break
    assert len(grid_lens)==len(tint['segs'])

    tint_s,tint_e = tint['segs'][0][0],tint['segs'][-1][-1]
    plot_tids = list()
    for tid,transcript in transcripts.items():
        if transcript['chrom']!=tint['chrom']:
            continue
        s,e = transcript['intervals'][0][0],transcript['intervals'][-1][-1]
        if not (tint_e<s or e<tint_s):
            plot_tids.append(tid)
    plot_tids.sort(key=lambda tid: transcripts.get(tid, {'name':tid})['name'])

    tid_colors = dict()
    color_idx = 0
    for tid in plot_tids:
        if not tid in tint['tids']:
            tid_colors[tid]='gray'
        else:
            tid_colors[tid]=colors[color_idx%len(colors)]
            color_idx+=1

    plot_settings['segs']=tint['segs']
    plot_settings['grid_lens']=grid_lens
    plot_settings['plot_tids']=plot_tids
    plot_settings['tid_colors']=tid_colors
    for isoform in tint['isoforms'].values():
        plot_isoform(
            isoform=isoform,
            transcripts=transcripts,
            plot_settings=plot_settings,
            outpath='{}/{}.pdf'.format(out_dir, isoform['id']),
        )


def main():
    args = parse_args()

    transcripts = get_transcripts(gtf=args.annotation_gtf)
    # for tid,transcript in transcripts.items():
    #     print(tid)
    #     print(transcript['intervals'])
    #     print(transcript['name'])
    plot_args = list()
    for tint in get_tints(cluster_tsv=args.cluster_tsv, segment_tsv=args.segment_tsv).values():
        plot_args.append((
            tint,
            transcripts,
            '{}/{}'.format(args.out_dir, tint['id']),
        ))
    #     print(tint['id'])
    #     print(tint['chrom'])
    #     print(tint['segs'])
    #     print(tint['tids'])
    #     print(tint.keys())
    #     for isoform in tint['isoforms'].values():
    #         print(isoform['id'])
    #         print(isoform['data'])
    #         print(isoform.keys())
    #         for read in isoform['reads']:
    #             print(read)

    os.makedirs(args.out_dir, exist_ok=True)
    if args.threads > 1:
        with Pool(args.threads) as p:
            for _ in p.imap_unordered(plot_tint, plot_args, chunksize=20):
                pass
    else:
        for _ in map(plot_tint, plot_args):
            pass

if __name__ == "__main__":
    main()
