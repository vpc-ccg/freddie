#!/usr/bin/env python3
from multiprocessing import Pool
import os
import argparse
import re
from operator import itemgetter
from itertools import groupby

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
    (0, 2)
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
    parser.add_argument("--tints",
                        type=int,
                        nargs='+',
                        default=list(),
                        help="Keep only specified tint ids. Default: All tints")
    parser.add_argument("-od",
                        "--out-dir",
                        type=str,
                        default='freddie_plot',
                        help="Output directory. Will be created if does not exist. Default: freddie_plot")
    args = parser.parse_args()
    args.out_dir = args.out_dir.rstrip('/')
    args.tints = set(args.tints)
    assert 1 <= args.threads <= 256
    return args

def plot_isoform(isoform, transcripts, plot_settings, title, outpath):
    grid_lens = plot_settings['grid_lens']
    segs      = plot_settings['segs']
    seg_idxs  = plot_settings['seg_idxs']
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
        height_ratios=[1,5],
        width_ratios=grid_lens,
        hspace=.1,
        wspace=0
    )
    t_axes = [fig.add_subplot(gs[0,i]) for i in range(len(grid_lens))]
    t_axes[len(t_axes)//2].set_title(title)
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
        ax.set_ylim(0,ylim)
        ax.set_yticks([])
        ax.set_xlim(0, 150)
    first_pos = segs[0][0]
    for axes,ylim in [(t_axes,len(plot_settings['plot_tids'])+1),(r_axes,len(isoform['reads']))]:
        for ax,seg_idx in zip(axes,seg_idxs):
            (s,e) = segs[seg_idx]
            ax.grid(zorder=0)
            ax.set_ylim(0,ylim)
            ax.set_xlim(0, 1)
            ax.set_xticks([0])
            ax.set_xticklabels([s-first_pos], rotation=45)
            ax.set_yticks(range(0,ylim,max(ylim//20,1)))
            if ax.is_first_col():
                ax.set_xticklabels([first_pos], rotation=45)
                pass
            elif ax.is_last_col():
                ax.yaxis.tick_right()
                ax.set_xticks([0,1])
                ax.set_xticklabels([s,e], rotation=45)
            else:
                ax.set_yticklabels([])
    t_axes[0].set_yticklabels( [transcripts[tid]['name'] for tid in plot_settings['plot_tids']]+['{}'.format('Consensus')])
    t_axes[-1].set_yticklabels([transcripts[tid]['name'] for tid in plot_settings['plot_tids']]+['{}'.format('Consensus')])

    for p,read in enumerate(sorted(isoform['reads'], key=itemgetter('tid', 'data'))):
        if read['tid'] in plot_settings['tid_colors']:
            color = plot_settings['tid_colors'][read['tid']]
        else:
            color = 'gray'
        for s,e in read['intervals']:
            x_0 = -1
            x_1 = -1
            for r_ax,seg_idx in zip(r_axes,seg_idxs):
                (seg_s,seg_e) = segs[seg_idx]
                if not (s <= seg_e and e >= seg_s):
                    continue
                seg_l = seg_e-seg_s
                if seg_s < s:
                    x_0 = (s-seg_s)/seg_l
                else:
                    x_0 = 0.0
                if seg_e < e:
                    x_1 = 1.0
                else:
                    x_1 = (e-seg_s)/seg_l
                r_ax.add_patch(patches.Rectangle(xy=(x_0,p),width=x_1,height=1,color=color))
        # for aid,ax in enumerate(r_axes):
        #     j = seg_idxs[aid]
        #     if read['gaps'][j] > 0:
        #         ax.text(0.9,p+0.5,s='{}'.format(read['gaps'][j]), size='xx-small', ha='right', color='red' if read['gaps'][j]>99 else 'black')
        #     if read['data'][j]=='0':
        #         pass
        #     if read['data'][j]=='1':
        #         ax.add_patch(patches.Rectangle(
        #             xy     = (0,p),
        #             width  = 1,
        #             height = 1,
        #             color  = color,
        #         ))
        #     if read['data'][j]=='2':
        #         ax.add_patch(patches.Rectangle(
        #             xy     = (0,p),
        #             width  = 1,
        #             height = 1,
        #             facecolor  = color,
        #             hatch     = '/',
        #         ))
    # for ax in r_axes+[sce_ax,scs_ax]:
    #     ax.hlines(y=ch_changes, xmin=0, xmax=1, alpha=0.4)
    for p,tid in enumerate(plot_settings['plot_tids']):
        for enum,(s,e) in zip(transcripts[tid]['enum'],transcripts[tid]['intervals']):
            x_0 = -1
            x_1 = -1
            temp = list()
            for t_ax,seg_idx in zip(t_axes,seg_idxs):
                (seg_s,seg_e) = segs[seg_idx]
                if not (s <= seg_e and e >= seg_s):
                    continue
                seg_l = seg_e-seg_s
                if seg_s < s:
                    x_0 = (s-seg_s)/seg_l
                else:
                    x_0 = 0.0
                if seg_e < e:
                    x_1 = 1.0
                else:
                    x_1 = (e-seg_s)/seg_l
                temp.append(t_ax)
                t_ax.add_patch(patches.Rectangle(xy=(x_0,p),width=x_1,height=1,color=plot_settings['tid_colors'][tid]))
            temp[len(temp)//2].text(x=0.5, y=p+0.5, s=enum, horizontalalignment='center', verticalalignment='center')
    for s,e in isoform['intervals']:
        x_0 = -1
        x_1 = -1
        temp = list()
        for t_ax,seg_idx in zip(t_axes,seg_idxs):
            (seg_s,seg_e) = segs[seg_idx]
            if not (s <= seg_e and e >= seg_s):
                continue
            seg_l = seg_e-seg_s
            if seg_s < s:
                x_0 = (s-seg_s)/seg_l
            else:
                x_0 = 0.0
            if seg_e < e:
                x_1 = 1.0
            else:
                x_1 = (e-seg_s)/seg_l
            temp.append(t_ax)
            t_ax.add_patch(patches.Rectangle(xy=(x_0,len(plot_settings['plot_tids'])),width=x_1,height=1,color='black'))
    p=len(plot_settings['plot_tids'])
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
                enum=list(),
                name=re.search(r'transcript_name \"(?P<tname>[\w\.-]+)\"', line[8]).group('tname')
            )
        transcripts[tid]['intervals'].append(interval)
        transcripts[tid]['enum'].append(re.search(r'exon_number \"(?P<enum>[\w\.-]+)\"', line[8]).group('enum'))
    return transcripts

def informative_segs(seg_content):
    M = len(seg_content)
    informative = [True for _ in range(M)]
    for j in range(1,M-1):
        if len(seg_content[j]) > 1:
            continue
        if '2' in seg_content[j]:
            continue
        if not (seg_content[j-1]==seg_content[j]==seg_content[j+1]):
            continue
        informative[j]=False
    return informative

def get_tints(cluster_tsv, segment_tsv, tint_ids=set()):
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
            if len(tint_ids)>0 and not tint_id in tint_ids:
                continue
            segs = segs.split(',')
            segs = [(int(s),int(e)) for s,e in zip(segs[:-1], segs[1:])]
            tints[tint_id]=dict(
                id=tint_id,
                chrom=chrom,
                segs=segs,
                partitions=dict()
            )
        elif line.startswith('isoform_'):
            continue
        else:
            line=line.rstrip().split('\t')
            rid=int(line[0])
            tint=int(line[4])
            if len(tint_ids)>0 and not tint_id in tint_ids:
                continue
            pid=int(line[5])
            iid = line[7]
            if iid =='*':
                iid = 'garbage'
            if not pid in tints[tint]['partitions']:
                tints[tint]['partitions'][pid] = dict(
                    id=pid,
                    tids=set(),
                    isoforms=dict(),
                    seg_idxs=list()
                )
            if not iid in tints[tint]['partitions'][pid]['isoforms']:
                tints[tint]['partitions'][pid]['isoforms'][iid] = dict(
                    id=iid,
                    reads=list(),
                )
            data=rid_to_data[rid]
            tid = line[1].split('_')[0]
            if len(tid) == 15 and tid.startswith('ENST') and tid[4:].isdigit():
                pass
            else:
                tid = None
            tints[tint]['partitions'][pid]['isoforms'][iid]['reads'].append(dict(
                rid=rid,
                name=line[1],
                tid=tid,
                chrom=line[2],
                strand=line[3],
                tint=tint,
                pid=pid,
                poly_tail_category=line[6],
                iid=iid,
                data=data,
                gaps=[int(x[:-1].split('(')[1]) if '(' in x else 0 for x in line[8+1:8+1+len(data)]],
                poly_tail=line[8+1+len(data):],
            ))
    for tint_id in tints.keys():
        for pid in tints[tint_id]['partitions'].keys():
            M = len(tints[tint_id]['segs'])
            seg_content = [set() for _ in range(M)]
            for isoform in tints[tint_id]['partitions'][pid]['isoforms'].values():
                cons = [0 for _ in range(M)]
                for read in isoform['reads']:
                    tints[tint_id]['partitions'][pid]['tids'].add(read['tid'])
                    read['intervals']=get_intervals(tints[tint_id]['segs'], read['data'])
                    for j in range(M):
                        seg_content[j].add(read['data'][j])
                        cons[j]+= (read['data'][j]=='1')
                isoform['cons']=['1' if x/len(isoform['reads'])>0.3 else '0' for x in cons]
                isoform['intervals']=get_intervals(tints[tint_id]['segs'], isoform['cons'])
            informative = informative_segs(seg_content)
            tints[tint_id]['partitions'][pid]['seg_idxs'] = [
                j for j in range(M) if informative[j]
            ]
    return tints

def get_intervals(segs,data):
    intervals=list()
    for d,group in groupby(enumerate(data), lambda x: x[1]):
        if d != '1':
            continue
        group = list(group)
        f_seg_idx = group[0][0]
        l_seg_idx = group[-1][0]
        intervals.append((segs[f_seg_idx][0],segs[l_seg_idx][1],))
    return intervals

def get_plot_tids(tint,partition,transcripts):
    plot_tids = list()
    tint_s,tint_e = tint['segs'][0][0],tint['segs'][-1][-1]
    for tid,transcript in transcripts.items():
        if not transcript['name'] in ['AR-UNION', 'AR-204', 'AR-201', 'AR-208']:
            continue
        if transcript['chrom']!=tint['chrom']:
            continue
        tid_s,tid_e = transcript['intervals'][0][0],transcript['intervals'][-1][-1]
        if tint_e<tid_s or tid_e<tint_s:
            continue
        found = False
        for int_s,int_e in transcripts[tid]['intervals']:
            for seg_idx in partition['seg_idxs']:
                (seg_s,seg_e) = tint['segs'][seg_idx]
                if seg_s <= int_s <= seg_e or seg_s <= int_e <= seg_e :
                    found = True
                if found:
                    break
            if found:
                break
        if found:
            plot_tids.append(tid)
    plot_tids.sort(key=lambda tid: transcripts.get(tid, {'name':tid})['name'])
    return plot_tids

def plot_partition(tint, partition, transcripts, out_dir):
    plot_settings = dict()

    partition['seg_idxs'] = list(range(len(tint['segs'])))
    plot_tids = get_plot_tids(tint,partition,transcripts)

    new_segs={tint['segs'][0][0],tint['segs'][-1][-1]}
    for tid in plot_tids:
        for s,e in transcripts[tid]['intervals']:
            new_segs.add(s)
            new_segs.add(e)
    new_segs=sorted(new_segs)

    tint['segs'] = [(s,e) for s,e in zip(new_segs[:-1],new_segs[1:])]
    partition['seg_idxs'] = list(range(len(tint['segs'])))

    grid_lens = list()
    for seg_idx in partition['seg_idxs']:
        s,e = tint['segs'][seg_idx]
        for threshold,value in grid_width_ratios:
            if e-s > threshold:
                grid_lens.append(value)
                break

    tid_colors = dict()
    color_idx = 0
    for tid in plot_tids:
        if tid in partition['tids'] or None in partition['tids']:
            tid_colors[tid]=colors[color_idx%len(colors)]
            color_idx+=1
        else:
            tid_colors[tid]='gray'

    plot_settings['segs']=tint['segs']
    plot_settings['seg_idxs']=partition['seg_idxs']
    plot_settings['grid_lens']=grid_lens
    plot_settings['plot_tids']=plot_tids
    plot_settings['tid_colors']=tid_colors
    for isoform in partition['isoforms'].values():
        plot_isoform(
            isoform=isoform,
            transcripts=transcripts,
            plot_settings=plot_settings,
            title='tint {}, partition {}, isoform {} (n={})'.format(tint['id'], partition['id'], isoform['id'], len(isoform['reads'])),
            outpath='{}/{}.{}.pdf'.format(out_dir, partition['id'], isoform['id']),
        )
    mergedObject = PdfFileMerger()
    for isoform in partition['isoforms'].values():
        mergedObject.append(PdfFileReader('{}/{}.{}.pdf'.format(out_dir, partition['id'], isoform['id']), 'rb'))
    mergedObject.write('{}/{}.pdf'.format(out_dir, partition['id'],))
    print('Merged all tint {} partition {} into {}/{}.pdf'.format(tint['id'], partition['id'], out_dir, partition['id'],))
    for isoform in partition['isoforms'].values():
        os.remove('{}/{}.{}.pdf'.format(out_dir, partition['id'], isoform['id']))

def plot_tint(plot_args):
    tint,transcripts,out_dir = plot_args
    print('Outputting to',out_dir)
    os.makedirs(out_dir, exist_ok=True)
    for partition in tint['partitions'].values():
        print('Plotting tint {} partition {}'.format(tint['id'], partition['id']))
        plot_partition(tint=tint, partition=partition, transcripts=transcripts, out_dir=out_dir)

def main():
    args = parse_args()

    transcripts = get_transcripts(gtf=args.annotation_gtf)
    plot_args = list()
    for tint in get_tints(cluster_tsv=args.cluster_tsv, segment_tsv=args.segment_tsv, tint_ids=args.tints).values():
        plot_args.append((
            tint,
            transcripts,
            '{}/{}'.format(args.out_dir, tint['id']),
        ))
    os.makedirs(args.out_dir, exist_ok=True)
    if args.threads > 1:
        with Pool(args.threads) as p:
            for _ in p.imap_unordered(plot_tint, plot_args, chunksize=5):
                pass
    else:
        for _ in map(plot_tint, plot_args):
            pass

if __name__ == "__main__":
    main()
