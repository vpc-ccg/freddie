#!/usr/bin/env python3
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from statistics import mean

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster aligned reads into isoforms")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to PAF file of read alignments")
    parser.add_argument("-t",
                        "--transcripts",
                        type=str,
                        required=True,
                        help="Path to TSV file of transcript intervals")
    parser.add_argument("-e",
                        "--exons",
                        type=str,
                        required=True,
                        help="Path to TSV file of canonical exon intervals")
    parser.add_argument("-i",
                        "--isoforms",
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

def read_paf(paf, range_len=15):
    is_first = True
    pos_to_rid = list()
    read_name_to_id = dict()
    rid_to_intervals = dict()
    for line in open(paf):
        line = line.rstrip().split('\t')
        if is_first:
            t_len = int(line[6])
            t_name = line[5]
            is_first = False
            pos_to_rid = [set() for _ in range(t_len)]
        if t_len != int(line[6]) or t_name != line[5]:
            print("Multiple targets detected in PAF file!", file=stderr)
            print(line, file=stderr)
            exit(-1)
        rid = line[0]
        if not rid in read_name_to_id:
            rid_to_intervals[rid] = list()
        if any('oc:c:1' in tag for tag in line[12:]):
            t_start = max(0, int(line[7]) - 1)
            t_end = int(line[8]) + 1
            q_start = max(0, int(line[2]) - 1)
            q_end = int(line[3]) + 1
            t_interval = (t_start, t_end)
            q_interval = (q_start, q_end)
            rid_to_intervals[rid].append((t_interval, q_interval))
    for intervals in rid_to_intervals.values():
        intervals.sort()
    for rid,intervals in rid_to_intervals.items():
        new_intervals = list()
        for idx,(t_interval,q_interval) in enumerate(intervals):
            if idx == 0:
                new_intervals.append((t_interval,q_interval))
                continue
            (t_interval_prv,q_interval_prv) = new_intervals[-1]
            if t_interval[0] - t_interval_prv[1] < range_len and q_interval_prv[0] - q_interval_prv[1] < range_len:
                new_intervals[-1] = (
                    (t_interval_prv[0],t_interval[1]),
                    (q_interval_prv[0],q_interval[1]),
                )
            else:
                new_intervals.append((t_interval,q_interval))
        rid_to_intervals[rid] = new_intervals
    for rid,intervals in rid_to_intervals.items():
        for (t_start, t_end),(_, _) in intervals:
            for i in range(t_start, t_end):
                pos_to_rid[i].add(rid)
    return pos_to_rid,rid_to_intervals

def plot_isoforms(exons, transcripts, pos_to_rid, rid_to_intervals, matrix, iid_to_isoform, iid_to_rids, out_prefix):
    L = len(iid_to_isoform)
    fig, axes = plt.subplots(L, 1, sharex='col', sharey='row', figsize=(30,8*(L)), squeeze=False)
    cmap_ins = plt.get_cmap('gnuplot2_r')
    norm_ins = mpl.colors.Normalize(vmin=10, vmax=800)
    fig.suptitle('L = {} N = {}'.format(L, len(set.union(*pos_to_rid))))
    for ax0,(iid,rids) in zip(axes,iid_to_rids.items()):
        print('Plotting isoform {}'.format(iid, len(iid_to_rids)))
        true_isoforms = dict()
        isoform = iid_to_isoform[iid]
        ax0 = ax0[0]
        N = len(rids)
        ax0.set_title('L = {} N = {}'.format(iid, N))
        coverage = [len(pos_rids&rids) for pos_rids in pos_to_rid]
        ax0.plot(range(len(coverage)), coverage, color='green', zorder=50)

        top    = N*0.95
        bottom = N*0.05
        step = (top-bottom)/(N+1)
        # plot the isoform's exons and introns
        for eid,in_iso in enumerate(isoform):
            if in_iso:
                ax0.plot(exons[eid], [N,N], color='black', alpha=0.8, marker='o')
            else:
                ax0.plot(exons[eid], [N,N], color='gray',  alpha=0.2,)
        # Plot the reads
        eid_to_mistakes = [0 for _ in exons]
        for read_count,rid in enumerate(rids):
            true_isoforms[rid.split('_')[0]] = true_isoforms.get(rid.split('_')[0], 0) + 1
            if read_count%50 == 0:
                print('Processing read {}/{}'.format(read_count, len(rids)))
            h = top-step*read_count
            # plot read mistakes
            for eid in range(len(exons)):
                if matrix[rid][eid] == 'X':
                    eid_to_mistakes[eid]+=1
                    ax0.scatter(x=mean(exons[eid]), y=h, s=0.5, marker='2', color='red')
            # Plot full read alignments
            for idx in range(len(rid_to_intervals[rid])):
                (cur_t_start,cur_t_end), (cur_q_start,cur_q_end) = rid_to_intervals[rid][idx]
                ax0.plot([cur_t_start,cur_t_end], [h,h], color='black', lw=0.25, zorder=0, alpha=0.5)
                if idx - 1 >= 0:
                    (_,_), (_,lst_q_end)   = rid_to_intervals[rid][idx-1]
                else:
                    lst_q_end   = 0
                dist_to_lst = cur_q_start - lst_q_end
                ax0.scatter(x=cur_t_start, y=h, s=0.25, color=cmap_ins(norm_ins(dist_to_lst)), zorder=5)
                if idx + 1 < len(rid_to_intervals[rid]):
                    (_,_), (nxt_q_start,_) = rid_to_intervals[rid][idx+1]
                else:
                    nxt_q_start = exons[-1][1]
                dist_to_nxt = nxt_q_start - cur_q_end
                ax0.scatter(x=cur_t_end,   y=h, s=0.25, color=cmap_ins(norm_ins(dist_to_nxt)), zorder=5)
        true_isoforms = sorted(true_isoforms.items(), key=lambda kv: kv[1], reverse=True)
        # Plot the isoform's exon lines
        for eid,exon in enumerate(exons):
            ax0.vlines(exon, ymin=0, ymax=N, color='gray', linestyle='solid', lw=1, alpha=0.5)
            ax0.text(x=mean(exon), y=N*1.015, s='{}'.format(eid_to_mistakes[eid]))
        step = 0.1
        for idx,(tid,cnt) in enumerate(true_isoforms):
            h = N*(0.9-idx*step)

            ax0.text(x=10, y=h, s='{} ({})'.format(tid, cnt))
            for s,e in transcripts[tid]:
                ax0.hlines(h, xmin=s, xmax=e, color='orange')
    plt.savefig('{}.pdf'.format(out_prefix))

def read_isoforms(isoforms):
    iid_to_rids = dict()
    iid_to_isoform = dict()
    matrix = dict()
    iid = -1
    for line in open(isoforms):
        line = line.rstrip().split('\t')
        if line[0][0] == '#':
            iid                 = line[0][1:]
            iid_to_rids[iid]    = set()
            iid_to_isoform[iid] = [x == '1' for x in line[1]]
            continue
        rid         = line[1]
        matrix[rid] = [i.split('(')[0] for i in line[2:]]
        iid_to_rids[iid].add(rid)
    return matrix,iid_to_rids,iid_to_isoform

def read_exons(exons_tsv):
    exons = list()
    for line in open(exons_tsv):
        i = line.rstrip().split('\t')[0:2]
        i = [int(x) for x in i]
        exons.append(i)
    return exons

def read_transcripts(transcripts_tsv):
    transcripts = dict()
    for line in open(transcripts_tsv):
        line = line.rstrip().split('\t')
        tid = line[0]
        transcripts[tid] = list()
        for i in line[3].split(','):
            s,e = i.split('-')
            transcripts[tid].append((int(s),int(e)))
    return transcripts

def main():
    args = parse_args()

    exons                              = read_exons(exons_tsv=args.exons)
    transcripts                        = read_transcripts(transcripts_tsv=args.transcripts)
    pos_to_rid,rid_to_intervals        = read_paf(paf=args.paf)
    matrix, iid_to_rids,iid_to_isoform = read_isoforms(isoforms=args.isoforms)

    plot_isoforms(exons=exons, transcripts=transcripts, pos_to_rid=pos_to_rid,rid_to_intervals=rid_to_intervals, iid_to_isoform=iid_to_isoform, matrix=matrix, iid_to_rids=iid_to_rids, out_prefix=args.out_prefix)

if __name__ == "__main__":
    main()
