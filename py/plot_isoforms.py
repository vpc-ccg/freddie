#!/usr/bin/env python3
import os
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from PyPDF2 import PdfFileMerger, PdfFileReader
from dateutil import parser

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
colors_ch = [
    '#999999',
    '#e41a1c',
    '#377eb8',
    '#4daf4a',
    '#984ea3',
    '#ff7f00',
    '#ffff33',
    '#a65628',
    '#f781bf',
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
    parser.add_argument("-q",
                        "--fastq",
                        type=str,
                        required=True,
                        help="Path to FASTQ file of reads")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to PAF file of reads")
    parser.add_argument("-i",
                        "--isoforms-tsv",
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

def get_softclip(paf):
    rname_to_sc = dict()
    for line in open(paf):
        line = line.rstrip().split('\t')
        name = line[0]
        if not name in rname_to_sc:
            rname_to_sc[name] = dict(
                length  = int(line[1]),
                strand  = line[4],
                q_start = int(line[2]),
                q_end   = int(line[3]),
            )
        rname_to_sc[name]['q_start'] = min(rname_to_sc[name]['q_start'], int(line[2]))
        rname_to_sc[name]['q_end']   = max(rname_to_sc[name]['q_end'], int(line[3]))
        if rname_to_sc[name]['strand'] == '+':
            rname_to_sc[name]['scs']=rname_to_sc[name]['q_start']
            rname_to_sc[name]['sce']=rname_to_sc[name]['length']-rname_to_sc[name]['q_end']-1
        else:
            rname_to_sc[name]['scs']=rname_to_sc[name]['length']-rname_to_sc[name]['q_end']-1
            rname_to_sc[name]['sce']=rname_to_sc[name]['q_start']
    return rname_to_sc

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

def read_fastq(fastq):
    rname_to_info = dict()
    for line in open(fastq):
        if line[0]!='@':
            continue
        line = line[1:].rstrip().split()
        rname = line[0]
        info=dict(
            runid=0,
            read=0,
            ch=0,
            start_time="1970-01-01T00:00:00Z",
        )
        for x in line[1:]:
            x = x.split('=')
            if x[0] in info:
                info[x[0]] = x[1]
        info['ch'] = int(info['ch'])
        info['start_time'] = parser.parse(info['start_time']).timestamp()
        rname_to_info[rname]=info
    return rname_to_info

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
        tname = rname.split('_')[0]
        if len(tname)==len('ENST00000396043') and tname[0:4]=='ENST' and tname[5:].isdigit():
            pass
        else:
            tname=''
        reads[rid]=dict(iid=iid, name=rname, data=d, tname=tname)
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

    rid_to_p = {rid:p for p,(tname,ch,data,rid) in enumerate(sorted(
        [(
            reads[rid]['tname'],
            reads[rid]['data'],
            reads[rid]['ch'],
            rid
        ) for rid in isoform['rids']]
    ))}
    for aid,ax in enumerate(r_axes):
        for rid in isoform['rids']:
            if reads[rid]['data'][aid]==False:
                continue
            ax.add_patch(patches.Rectangle(
                xy     = (0,rid_to_p[rid]),
                width  = 1,
                height = 1,
                color  = tid_to_color.get(reads[rid]['tname'], colors_ch[reads[rid]['ch']%len(colors_ch)])
            ))

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
    rname_to_info = read_fastq(args.fastq)
    rname_to_sc   = get_softclip(args.paf)
    segs,grid_lens,tid_to_segs,tid_to_color=get_tinfo(
        segments_txt=args.segments_txt,
        transcripts_tsv=args.transcripts_tsv
    )
    isoforms,reads=get_isoforms(
        isoforms_tsv=args.isoforms_tsv
    )
    for rid,read in reads.items():
        rname = read['name']
        assert rname in rname_to_info
        for k,v in rname_to_info[rname].items():
            assert not k in read
            reads[rid][k]=v
        assert rname in rname_to_sc
        reads[rid][k]['scs'] = rname_to_sc[rname]['scs']
        reads[rid][k]['sce'] = rname_to_sc[rname]['sce']
        reads[rid][k]['strand'] = rname_to_sc[rname]['strand']
    # for read in reads.values():
    #     print(read)
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
