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
    parser.add_argument("-e",
                        "--exons",
                        type=str,
                        required=True,
                        help="Path to TSV file of canonical exon intervals")
    parser.add_argument("-d",
                        "--data-matrix",
                        type=str,
                        required=True,
                        help="Path to TXT file of reads binary vectors")
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
        name = line[0]
        if not name in read_name_to_id:
            rid = len(read_name_to_id)
            read_name_to_id[name] = rid
            rid_to_intervals[rid] = list()
        rid = read_name_to_id[name]
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

def plot_isoforms(exons, pos_to_rid, rid_to_intervals, isoform_to_rids, matrix, isoforms, out_prefix):
    L = len(isoform_to_rids)
    fig, axes = plt.subplots(L, 1, sharex='col', sharey='row', figsize=(30,8*(L)), squeeze=False)
    cmap_ins = plt.get_cmap('gnuplot2_r')
    norm_ins = mpl.colors.Normalize(vmin=10, vmax=800)
    fig.suptitle('L = {} N = {}'.format(L, len(set.union(*pos_to_rid))))
    for (isoform_id,rids),isoform in zip(isoform_to_rids.items(),isoforms):
        ax0 = axes[isoform_id][0]
        N = len(rids)
        ax0.set_title('N = {}'.format(N))
        coverage = [len(pos_rids&rids) for pos_rids in pos_to_rid]
        ax0.plot(range(len(coverage)), coverage, color='green', zorder=50)

        top    = N*0.95
        bottom = N*0.05
        step = (top-bottom)/N
        # plot the isoform's exons and introns
        for eid,in_iso in enumerate(isoform):
            print(eid,in_iso, len(isoform))
            print(len(exons))
            if in_iso:
                ax0.plot(exons[eid], [N,N], color='black', alpha=0.8, marker='o')
            else:
                ax0.plot(exons[eid], [N,N], color='gray',  alpha=0.2,)
        # Plot the reads
        eid_to_mistakes = [0 for _ in exons]
        for read_count,rid in enumerate(rids):
            h = top-step*read_count
            # plot read mistakes
            for eid,in_iso in enumerate(isoform):
                if in_iso and matrix[rid][eid] == 0:
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
        # Plot the isoform's exon lines
        for eid,exon in enumerate(exons):
            ax0.vlines(exon, ymin=0, ymax=N, color='gray', linestyle='solid', lw=1, alpha=0.5)
            ax0.text(x=mean(exon), y=N*1.015, s='{}'.format(eid_to_mistakes[eid]))
    plt.savefig('{}.pdf'.format(out_prefix))

def read_isoforms(isoforms):
    isoform_to_rids = dict()
    for line in open(isoforms):
        line = line.rstrip().split('\t')
        rid        = int(line[0])
        isoform_id = int(line[1])

        if not isoform_id in isoform_to_rids:
            isoform_to_rids[isoform_id] = set()
        isoform_to_rids[isoform_id].add(rid)
    return isoform_to_rids

def read_matrix(matrix_txt):
    rids_to_exons = list()
    for line in open(matrix_txt):
        line = line.rstrip()
        rids_to_exons.append([x == '1' for x in line])
        print(len(rids_to_exons[-1]))
    return rids_to_exons

def get_isoforms(isoform_to_rids, matrix):
    M = len(matrix[0])
    isoforms = [[False for _ in range(M)] for _ in range(len(isoform_to_rids))]
    for isoform_id,rids in isoform_to_rids.items():
        print(rids)
        for rid in rids:
            for eid in range(M):
                isoforms[isoform_id][eid] |=  matrix[rid][eid]
    return isoforms

def read_exons(exons_tsv):
    exons = list()
    for line in open(exons_tsv):
        i = line.rstrip().split('\t')[0:2]
        i = [int(x) for x in i]
        exons.append(i)
    return exons

def main():
    args = parse_args()

    exons = read_exons(exons_tsv=args.exons)
    print('Exons len:', len(exons))
    matrix = read_matrix(matrix_txt=args.data_matrix)
    pos_to_rid,rid_to_intervals = read_paf(paf=args.paf)
    isoform_to_rids = read_isoforms(isoforms=args.isoforms)
    isoforms = get_isoforms(isoform_to_rids=isoform_to_rids, matrix=matrix)
    plot_isoforms(exons=exons, pos_to_rid=pos_to_rid, matrix=matrix,rid_to_intervals=rid_to_intervals, isoforms=isoforms, isoform_to_rids=isoform_to_rids, out_prefix=args.out_prefix)

if __name__ == "__main__":
    main()
