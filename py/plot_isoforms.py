#!/usr/bin/env python3
import os
import argparse
from operator import itemgetter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from PyPDF2 import PdfFileMerger, PdfFileReader
import dateutil

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
    'gray',
    # '#fee0d2',
    # '#fcbba1',
    # '#fc9272',
    '#fb6a4a',
    '#de2d26',
    '#a50f15',
]
#fb6a4a
#de2d26
#a50f15
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
        tid,_,strand,intervals,_ = l.rstrip().split('\t')
        tid_to_segs[tid] = sorted([(int(x.split('-')[0]),int(x.split('-')[1])) for x in intervals.split(',')])
    tid_to_color = {tid:colors[tidx%len(colors)] for tidx,tid in enumerate(sorted(tid_to_segs.keys()))}
    return segs,grid_lens,tid_to_segs,tid_to_color,strand

def read_fastq(fastq):
    rname_to_info = dict()
    for line_num,line in enumerate(open(fastq)):
        if line_num%4==1:
            rname_to_info[rname]['length']=len(line.rstrip())
        if line[0]!='@':
            continue
        line = line[1:].rstrip().split()
        rname = line[0]
        info=dict(
            runid=0,
            read=0,
            ch=-1,
            start_time="1970-01-01T00:00:00Z",
        )
        for x in line[1:]:
            x = x.split('=')
            if x[0] in info:
                info[x[0]] = x[1]
        info['ch'] = int(info['ch'])
        info['start_time'] = dateutil.parser.parse(info['start_time']).timestamp()
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
        if len(tname)==len('ENST00000396043') and tname.startswith('ENST') and tname[5:].isdigit():
            pass
        else:
            tname=''
        read=dict(iid=iid, name=rname, data=d, tname=tname)
        read['tail']=dict(
            S_poly_size  = 0,
            S_poly_gap  = 0,
            S_poly_type  = '',
            SSC = 0,
            E_poly_size  = 0,
            E_poly_gap  = 0,
            E_poly_type  = '',
            ESC = 0,
        )
        for gap in line[4+len(line[3]):]:
            k,v=gap.split(':')
            v = [int(x) for x in v.strip('()').split(',')]
            if v[0]<1:
                continue
            for s in ['S','E']:
                if gap.startswith('{}SC'.format(s)):
                    read['tail']['{}SC'.format(s)]=v[0]
                    continue
                for c in ['A','T']:
                    if gap.startswith('{}{}'.format(s,c)):
                        read['tail']['{}_poly_type'.format(s)]=c
                        read['tail']['{}_poly_size'.format(s)]=v[0]
                        read['tail']['{}_poly_gap'.format(s)]=v[1]
        assert len(read['tail'])==8,read['tail']
        reads[rid]=read
        isoforms[iid]['rids'].append(rid)
    isoforms = [isoforms[x] for x in sorted(isoforms.keys())]
    return isoforms,reads

def plot_isoforms(isoform, grid_lens, tid_to_segs, segs, reads, tid_to_color, strand, out_prefix):
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

    ylim=len(isoform['rids'])
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
    for axes,ylim in [(t_axes,len(tid_to_segs)),(r_axes,len(isoform['rids']))]:
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
    t_axes[0].set_yticklabels(sorted(tid_to_segs.keys()))
    t_axes[-1].set_yticklabels(sorted(tid_to_segs.keys()))
    isoform_reads = [reads[rid] for rid in isoform['rids']]
    # tid_to_color.get(read['tname'],'gray')
    for read in isoform_reads:
        read['color_idx'] = 0
    ch_to_rids = {read['ch']:list() for read in isoform_reads}
    for p,read in enumerate(isoform_reads):
        ch_to_rids[read['ch']].append(p)
    for ch,rids in ch_to_rids.items():
        if ch==-1:
            continue
        if len(rids)!=2:
            continue
        # print('ch size', len(rids))
        for r1 in rids:
            for r2 in rids:
                # if r2<=r1:
                #     continue
                if not any(isoform_reads[r1]['data']):
                    continue
                if not any(isoform_reads[r2]['data']):
                    continue
                # if not any(d1 and d2 for (d1,d2) in zip(isoform_reads[r1]['data'],isoform_reads[r2]['data'])):
                # if isoform_reads[r1]['data'].index(True) > isoform_reads[r2]['data'][::-1].index(True):
                #     isoform_reads[r1]['color_idx']+=1
                #     isoform_reads[r2]['color_idx']+=1
                #     print(ch,r1,r2,isoform_reads[r1]['data'].index(True),isoform_reads[r2]['data'][::-1].index(True))
    # for read in isoform_reads:
    #     ch_count[read['ch']]+=1
    # for read in isoform_reads:
    #     read['ch_order'] = ch_count[read['ch']]
    isoform_reads.sort(key=itemgetter(
        'tname',
        'data',
        # 'ch_order',
        # 'ch',
    ))
    # rid_to_p = {rid:p for p,(tname,ch,data,rid) in enumerate(sorted(
    #     [(
    #         reads[rid]['tname'],
    #         reads[rid]['data'],
    #         reads[rid]['ch'],
    #         rid
    #     ) for rid in isoform['rids']]
    # ))}
    # last_ch = None
    # ch_color_idx = -1
    # ch_subcolor_idx = -1
    # ch_changes = list()
    for p,read in enumerate(isoform_reads):
        # if ch_count[read['ch']] > 1:
        #     if last_ch != read['ch']:
        #         ch_changes.append(p)
        #         last_ch = read['ch']
        #         ch_color_idx = (ch_color_idx+1)%len(colors_ch)
        #         ch_subcolor_idx=-1
        #     ch_subcolor_idx = (ch_subcolor_idx+1)%len(colors_ch[ch_color_idx])
        #     color = colors_ch[ch_color_idx][ch_subcolor_idx]
        # else:
        if read['tname'] in tid_to_color:
            color = tid_to_color[read['tname']]
        else:
            color = 'gray'
            # if read['color_idx']<len(colors_ch):
            #     color = colors_ch[read['color_idx']]
            # else:
            #     color = colors_ch[-1]
        for s_ax,x in zip(s_axes,s_tail_order):
            k = 'S{}'.format(x)
            if x =='_poly_size' and (strand=='+' or read['tail']['S_poly_type']=='A'):
                s_ax.add_patch(patches.Rectangle(
                    xy     = (0,p),
                    width  = read['tail'][k],
                    height = 1,
                    color  = 'red',
                    hatch  = '/',
                ))
            else:
                s_ax.add_patch(patches.Rectangle(
                    xy     = (0,p),
                    width  = read['tail'][k],
                    height = 1,
                    color  = color,
                ))
        for e_ax,x in zip(e_axes,e_tail_order):
            k = 'E{}'.format(x)
            if x =='_poly_size' and (strand=='-' or read['tail']['E_poly_type']=='T'):
                e_ax.add_patch(patches.Rectangle(
                    xy     = (0,p),
                    width  = read['tail'][k],
                    height = 1,
                    color  = 'red',
                    hatch  = '/',
                ))
            else:
                e_ax.add_patch(patches.Rectangle(
                    xy     = (0,p),
                    width  = read['tail'][k],
                    height = 1,
                    color  = color,
                ))
        for aid,ax in enumerate(r_axes):
            if read['data'][aid]==False:
                continue
            ax.add_patch(patches.Rectangle(
                xy     = (0,p),
                width  = 1,
                height = 1,
                color  = color,
            ))
    # for ax in r_axes+[sce_ax,scs_ax]:
    #     ax.hlines(y=ch_changes, xmin=0, xmax=1, alpha=0.4)
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
            # t_axes[0].text(x=-0.5, y=tid_to_p[tid]+0.5, s=tid)
    # plt.savefig('{}.svg'.format(out_prefix),bbox_inches='tight')
    plt.savefig('{}.pdf'.format(out_prefix),bbox_inches='tight')

def main():
    args = parse_args()
    segs,grid_lens,tid_to_segs,tid_to_color,strand=get_tinfo(
        segments_txt=args.segments_txt,
        transcripts_tsv=args.transcripts_tsv
    )
    rname_to_info = read_fastq(args.fastq)
    if len(rname_to_info)==0:
        print('No reads for this gene!')
        out_file=open('{}.pdf'.format(args.out_prefix),'w+')
        out_file.close()
        exit()
    isoforms,reads=get_isoforms(
        isoforms_tsv=args.isoforms_tsv
    )
    # ch_to_rids = dict()
    for rid,read in reads.items():
        rname = read['name']
        assert rname in rname_to_info
        for k,v in rname_to_info[rname].items():
            assert not k in read
            reads[rid][k]=v
        read['rid']=rid
        try:
            read['fexon'] = read['data'].index(True)
            read['lexon']  = len(read['data'])-1-read['data'][::-1].index(True)
        except ValueError:
            read['fexon']=-1
            read['lexon']=-1
    #
    # time_sorted_reads = sorted((read for read in reads.values() if any(read['data'])), key=itemgetter('start_time'))
    # last_time = time_sorted_reads[0]['start_time']
    # last_fexon = time_sorted_reads[0]['fexon']
    # slack_time = time_sorted_reads[0]['length']/50
    # print('ch','strand','fexon', 'lexon', 'delta')
    # for read in time_sorted_reads:
    #     delta = read['start_time']-last_time
    #     if delta-slack_time < 1000 and last_fexon > read['lexon']:
    #         print(read['ch'],read['strand'],read['fexon'], read['lexon'], delta-slack_time, '<------')
    #     else:
    #         print(read['ch'],read['strand'],read['fexon'], read['lexon'], delta-slack_time)
    #     last_time  = read['start_time']
    #     last_fexon = read['fexon']
    #     slack_time = 0#read['length']
    # exit()
    for iid,isoform in enumerate(isoforms):
        print('Plotting isoform {}...'.format(iid))
        plot_isoforms(
            isoform=isoform,
            grid_lens=grid_lens,
            tid_to_segs=tid_to_segs,
            segs=segs,
            reads=reads,
            tid_to_color=tid_to_color,
            strand=strand,
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
