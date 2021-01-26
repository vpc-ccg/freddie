#!/usr/bin/env python3
import argparse
from itertools import groupby
import os
import glob


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract alignment information from BAM/SAM file and splits reads into distinct transcriptional intervals")
    parser.add_argument("-s",
                        "--segment-dir",
                        type=str,
                        required=True,
                        help="Path to directory of Freddie segment")
    parser.add_argument("-c",
                        "--cluster-dir",
                        type=str,
                        required=True,
                        help="Path to directory of Freddie cluster")
    parser.add_argument("-t",
                        "--seqpare-threshold",
                        type=float,
                        default=0.95,
                        help="Seqpare threshold to merge two isoforms")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default='freddie_isoforms.gtf',
                        help="Path to output file. Default: freddie_isoforms.gtf")
    args = parser.parse_args()
    assert 0.0 < args.seqpare_threshold <= 1.0
    return args


def get_tints(cluster_dir, segment_dir):
    tints = dict()
    rid_to_data = dict()
    for contig in os.listdir(segment_dir):
        if not os.path.isdir('{}/{}'.format(segment_dir, contig)):
            continue
        rid_to_data[contig] = dict()
        for segment_tsv in glob.iglob('{}/{}/segment_*.tsv'.format(segment_dir, contig)):
            for line in open(segment_tsv):
                if line[0] == '#':
                    continue
                line = line.rstrip().split('\t')
                rid_to_data[contig][int(line[0])] = line[5]
    for contig in os.listdir(cluster_dir):
        if not os.path.isdir('{}/{}'.format(cluster_dir, contig)):
            continue
        tints[contig] = dict()
        for cluster_tsv in glob.iglob('{}/{}/cluster_*.tsv'.format(cluster_dir, contig)):
            for line in open(cluster_tsv):
                if line.startswith('#'):
                    chrom, tint_id, segs = line.rstrip()[1:].split('\t')
                    assert contig == chrom
                    tint_id = int(tint_id)
                    segs = segs.split(',')
                    segs = [(int(s), int(e))
                            for s, e in zip(segs[:-1], segs[1:])]
                    tints[contig][tint_id] = dict(
                        id=tint_id,
                        chrom=chrom,
                        segs=segs,
                        partitions=dict()
                    )
                elif line.startswith('isoform_'):
                    continue
                else:
                    line = line.rstrip().split('\t')
                    rid = int(line[0])
                    tint = int(line[4])
                    pid = int(line[5])
                    iid = line[7]
                    if iid == '*':
                        iid = 'garbage'
                    if not pid in tints[contig][tint]['partitions']:
                        tints[contig][tint]['partitions'][pid] = dict(
                            id=pid,
                            tids=set(),
                            isoforms=dict(),
                            seg_idxs=list()
                        )
                    if not iid in tints[contig][tint]['partitions'][pid]['isoforms']:
                        tints[contig][tint]['partitions'][pid]['isoforms'][iid] = dict(
                            id=iid,
                            reads=list(),
                        )
                    data = rid_to_data[contig][rid]
                    tints[contig][tint]['partitions'][pid]['isoforms'][iid]['reads'].append(dict(
                        rid=rid,
                        name=line[1],
                        tid=line[1].split('_')[0],
                        chrom=line[2],
                        strand=line[3],
                        tint=tint,
                        pid=pid,
                        poly_tail_category=line[6],
                        iid=iid,
                        data=data,
                        gaps=[int(
                            x[:-1].split('(')[1]) if '(' in x else 0 for x in line[8+1:8+1+len(data)]],
                        poly_tail=line[8+1+len(data):],
                    ))
    for contig in tints.keys():
        for tint_id in tints[contig].keys():
            for pid in tints[contig][tint_id]['partitions'].keys():
                M = len(tints[contig][tint_id]['segs'])
                seg_content = [set() for _ in range(M)]
                for isoform in tints[contig][tint_id]['partitions'][pid]['isoforms'].values():
                    for read in isoform['reads']:
                        tints[contig][tint_id]['partitions'][pid]['tids'].add(read['tid'])
                        for j in range(M):
                            seg_content[j].add(read['data'][j])
    return tints


def gtf_intervals(tints):
    for contig in tints.keys():
        for tint in tints[contig].values():
            for partition in tint['partitions'].values():
                for isoform in partition['isoforms'].values():
                    isoform['cons'] = [False for _ in tint['segs']]
                    isoform['strand'] = '.'
                    isoform['intervals'] = list()
                    if isoform['id'] == 'garbage':
                        continue
                    cons = [0 for _ in tint['segs']]
                    poly_tail_categories = {'N': 0, 'S': 0, 'E': 0}
                    for read in isoform['reads']:
                        for j in range(len(tint['segs'])):
                            cons[j] += read['data'][j] == '1'
                            poly_tail_categories[read['poly_tail_category']] += 1
                    cons = [x/len(isoform['reads']) > 0.3 for x in cons]
                    if not True in cons:
                        continue
                    isoform['cons'] = cons
                    if poly_tail_categories['S'] > poly_tail_categories['E']:
                        isoform['strand'] = '-'
                    else:
                        isoform['strand'] = '+'
                    for d, group in groupby(enumerate(cons), lambda x: x[1]):
                        if d != True:
                            continue
                        group = list(group)
                        f_seg_idx = group[0][0]
                        l_seg_idx = group[-1][0]
                        isoform['intervals'].append((
                            tint['segs'][f_seg_idx][0]+1,
                            tint['segs'][l_seg_idx][1]
                        ))


def overlap(a, b):
    e = min(a[1], b[1])
    s = max(a[0], b[0])
    return max(0, e-s)


def seqpare(A, B):
    if (not 'intervals' in A) or (not 'intervals' in B):
        return 0.0
    scores = list()
    for i, a in enumerate(A['intervals']):
        for j, b in enumerate(B['intervals']):
            o = overlap(a, b)
            la = a[1]-a[0]
            lb = b[1]-b[0]
            s = o/(la+lb-o)
            scores.append((s, i, j))
    scores.sort(reverse=True)
    i_selected = set()
    j_selected = set()
    O = 0.0
    for s, i, j in scores:
        if i in i_selected or j in j_selected:
            continue
        O += s
        i_selected.add(i)
        j_selected.add(j)

    return O/(len(A['intervals'])+len(B['intervals'])-O)


def seqpare_matrix(tints):
    for tint in tints.values():
        for partition in tint['partitions'].values():
            partition['scores'] = dict()
            for a, A in partition['isoforms'].items():
                for b, B in partition['isoforms'].items():
                    if a >= b:
                        continue
                    s = seqpare(A, B)
                    partition['scores'][(a, b)] = s


def output_gtf(tints, outpath):
    out_file = open(outpath, 'w+')
    out_reads = open(outpath+'.reads.tsv', 'w+')
    for contig in tints.keys():
        for tint in tints[contig].values():
            for partition in tint['partitions'].values():
                for isoform in partition['isoforms'].values():
                    if len(isoform['intervals']) == 0:
                        continue
                    transcript_name = '{chr}_{tint}_{iid}'.format(
                        chr=tint['chrom'],
                        tint=tint['id'],
                        iid=isoform['id'],
                    )
                    for r in isoform['reads']:
                        print('{}\t{}'.format(
                            r['name'], transcript_name), file=out_reads)
                    record = list()
                    record.append(tint['chrom'])
                    record.append('freddie')
                    record.append('transcript')
                    record.append(str(isoform['intervals'][0][0]))
                    record.append(str(isoform['intervals'][-1][1]))
                    record.append('.')
                    record.append(isoform['strand'])
                    record.append('.')
                    record.append('transcript_id "{transcript_name}"; read_support "{read_support}";'.format(
                        transcript_name=transcript_name,
                        read_support=len(isoform['reads']),
                    ))
                    out_file.write('\t'.join(record))
                    out_file.write('\n')
                    for eid, (s, e) in enumerate(isoform['intervals'], start=1):
                        record = list()
                        record.append(tint['chrom'])
                        record.append('freddie')
                        record.append('exon')
                        record.append(str(s))
                        record.append(str(e))
                        record.append('.')
                        record.append(isoform['strand'])
                        record.append('.')
                        record.append('transcript_id "{transcript_name}"; exon_number "{eid}"; exon_id "{transcript_name}_{eid}"; '.format(
                            transcript_name=transcript_name,
                            eid=eid,
                        ))
                        out_file.write('\t'.join(record))
                        out_file.write('\n')
                        eid += 1
    out_file.close()
    out_reads.close()


def connected_components(matrix, t):
    # print('======')
    iid_to_cid = dict()
    comps = dict()
    for (a, b), s in matrix.items():
        comps[a] = [a]
        iid_to_cid[a] = a
        comps[b] = [b]
        iid_to_cid[b] = b
    for (a, b), s in matrix.items():
        if not s > t:
            continue
        a_cid = iid_to_cid[a]
        b_cid = iid_to_cid[b]
        if a_cid == b_cid:
            continue
        # print('a', a, end=' ')
        # print(comps[a_cid],)
        # print('b', b, end=' ')
        # print(comps[b_cid],)
        for x in comps[b_cid]:
            # print('{}: {} -> {}'.format(x, iid_to_cid[x], a_cid))
            iid_to_cid[x] = a_cid
        comps[a_cid].extend(comps.pop(b_cid))
    return comps


def merge_isoforms(tints, t):
    for tint in tints.values():
        for partition in tint['partitions'].values():
            if len(partition['scores']) == 0:
                continue
            cc = connected_components(partition['scores'], t)
            # print('pid', partition['id'], cc)
            for main_iid, other_iids in cc.items():
                if len(other_iids) == 1:
                    continue
                # print(main_iid)
                for iid in other_iids:
                    if iid == main_iid:
                        continue
                    other_reads = partition['isoforms'][iid]['reads']
                    for read in other_reads:
                        read['iid'] = main_iid
                    partition['isoforms'][main_iid]['reads'].extend(
                        other_reads)
                    partition['isoforms'].pop(iid)
    gtf_intervals(tints)


def main():
    args = parse_args()

    tints = get_tints(cluster_dir=args.cluster_dir,
                      segment_dir=args.segment_dir)
    gtf_intervals(tints)
    # seqpare_matrix(tints)
    # merge_isoforms(tints, t=args.seqpare_threshold)
    output_gtf(tints, args.output)


if __name__ == "__main__":
    main()
