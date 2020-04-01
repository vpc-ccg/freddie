#!/usr/bin/env python3
import os
import argparse
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
    parser.add_argument("-t",
                        "--transcripts-tsv",
                        type=str,
                        required=True,
                        help="Path to TSV file of transcript intervals")
    parser.add_argument("-s",
                        "--segments-txt",
                        type=str,
                        required=True,
                        help="Path to TXT file of canonical exons")
    parser.add_argument("-i",
                        "--isoforms_tsv",
                        type=str,
                        required=True,
                        help="Path to TSV file of read isoform assignments")
    parser.add_argument("-op",
                        "--out-prefix",
                        type=str,
                        required=True,
                        help="Output prefix that does not include .EXT part")
    args = parser.parse_args()
    return args

def get_tinfo(segments_txt, transcripts_tsv):
    segs = [int(l.strip()) for l in open(segments_txt)]
    segs = list(zip(segs[:-1],segs[1:]))
    grid_lens = list()
    for s,e in segs:
        for threshold,value in grid_width_ratios:
            if e-s > threshold:
                grid_lens.append(value)
                break
    assert len(grid_lens)==len(segs)
    tid_to_segs = dict()
    for l in open(transcripts_tsv):
        tid,_,_,intervals,_ = l.rstrip().split('\t')
        tid_to_segs[tid] = sorted([(int(x.split('-')[0]),int(x.split('-')[1])) for x in intervals.split(',')])
    tid_to_color = {tid:colors[tidx%len(colors)] for tidx,tid in enumerate(sorted(tid_to_segs.keys()))}
    return segs,grid_lens,tid_to_segs,tid_to_color

def get_isoforms(isoforms_tsv):
    reads = dict()
    isoforms = dict()
    iid = -1
    for line in open(isoforms_tsv):
        line = line.rstrip().split('\t')
        if line[0][0]=='#':
            iid+=1
            isoform_name = line[0][1:]
            isoform_data = [x=='1' for x in line[3]]
            isoforms[iid] = dict(name=isoform_name, data=isoform_data, rids=list())
            continue
        rname = line[1]
        rid = int(line[2])
        d = [x=='1' for x in line[3]]
        reads[rid]=dict(iid=iid, name=rname, data=d, tname=rname.split('_')[0])
        isoforms[iid]['rids'].append(rid)
    isoforms = [isoforms[x] for x in sorted(isoforms.keys())]
    return isoforms,reads

def plot_isoforms(isoform, grid_lens, tid_to_segs, segs, reads, tid_to_color, out_prefix):
    fig = plt.figure(figsize=(len(grid_lens), 30), constrained_layout=False)
    gs = gridspec.GridSpec(ncols=len(grid_lens), nrows=2, figure=fig, width_ratios=grid_lens, height_ratios=[1,5], hspace=.1, wspace=0)
    t_axes = [fig.add_subplot(gs[0,i]) for i in range(len(grid_lens))]
    r_axes = [fig.add_subplot(gs[1,i]) for i in range(len(grid_lens))]

    for axes,ylim in [(t_axes,len(tid_to_segs)),(r_axes,len(isoform['rids']))]:
        for ax,(s,e) in zip(axes,segs):
            ax.set_ylim(0,ylim)
            ax.set_xlim(0, 1)
            ax.set_xticks([0])
            ax.set_xticklabels([s], rotation=45)
            if ax.is_first_col():
                ax.set_yticks(range(0,ylim,max(ylim//20,1)))
            elif ax.is_last_col():
                ax.yaxis.tick_right()
                ax.set_yticks(range(0,ylim,max(ylim//20,1)))
                ax.set_xticks([0,1])
                ax.set_xticklabels([s,e], rotation=45)
            else:
                ax.set_yticks([])

    rid_to_p = {rid:p for p,(_,rid) in enumerate(sorted([(reads[rid]['name'],rid) for rid in isoform['rids']]))}
    for aid,ax in enumerate(r_axes):
        for rid in isoform['rids']:
            if reads[rid]['data'][aid]==False:
                continue
            ax.add_patch(patches.Rectangle((0,rid_to_p[rid]),1,1,color=tid_to_color.get(reads[rid]['tname'],'gray')))

    tid_to_p = {tid:p for p,tid in enumerate(sorted(tid_to_segs.keys()))}
    for tid,intervals in tid_to_segs.items():
        for s,e in intervals:
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
                t_axes[t_aid].add_patch(patches.Rectangle(xy=(x_0,tid_to_p[tid]),width=x_1,height=1,color=tid_to_color.get(tid,'gray')))
    # plt.savefig('{}.svg'.format(out_prefix),bbox_inches='tight')
    plt.savefig('{}.pdf'.format(out_prefix),bbox_inches='tight')

def main():
    args = parse_args()
    segs,grid_lens,tid_to_segs,tid_to_color=get_tinfo(
        segments_txt=args.segments_txt,
        transcripts_tsv=args.transcripts_tsv
    )
    isoforms,reads=get_isoforms(
        isoforms_tsv=args.isoforms_tsv
    )

    for iid,isoform in enumerate(isoforms):
        print('Plotting isoform {}...'.format(iid))
        plot_isoforms(
            isoform=isoform,
            grid_lens=grid_lens,
            tid_to_segs=tid_to_segs,
            segs=segs,
            reads=reads,
            tid_to_color=tid_to_color,
            out_prefix='{}.{}'.format(args.out_prefix,iid),
        )

    mergedObject = PdfFileMerger()
    for iid in range(len(isoforms)):
        mergedObject.append(PdfFileReader('{}.{}.pdf'.format(args.out_prefix,iid), 'rb'))
    mergedObject.write('{}.pdf'.format(args.out_prefix))
    for iid in range(len(isoforms)):
        os.remove('{}.{}.pdf'.format(args.out_prefix,iid))

if __name__ == "__main__":
    main()
