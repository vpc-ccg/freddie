#!/usr/bin/env python3
import argparse
from itertools import groupby

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract alignment information from BAM/SAM file and splits reads into distinct transcriptional intervals")
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
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default='freddie_isoforms.gtf',
                        help="Path to output file. Default: freddie_isoforms.gtf")
    args = parser.parse_args()
    return args

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
                partitions=dict()
            )
        elif line.startswith('isoform_'):
            continue
        else:
            line=line.rstrip().split('\t')
            rid=int(line[0])
            tint=int(line[4])
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
            tints[tint]['partitions'][pid]['isoforms'][iid]['reads'].append(dict(
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
                gaps=[int(x[:-1].split('(')[1]) if '(' in x else 0 for x in line[8+1:8+1+len(data)]],
                poly_tail=line[8+1+len(data):],
            ))
    for tint_id in tints.keys():
        for pid in tints[tint_id]['partitions'].keys():
            M = len(tints[tint_id]['segs'])
            seg_content = [set() for _ in range(M)]
            for isoform in tints[tint_id]['partitions'][pid]['isoforms'].values():
                for read in isoform['reads']:
                    tints[tint_id]['partitions'][pid]['tids'].add(read['tid'])
                    for j in range(M):
                        seg_content[j].add(read['data'][j])
    return tints

def output_gtf(tints, outpath):
    out_file=open(outpath,'w+')
    for tint in tints.values():
        for partition in tint['partitions'].values():
            for isoform in partition['isoforms'].values():
                if isoform['id']=='garbage':
                    continue
                cons = [0 for _ in tint['segs']]
                for j in range(len(tint['segs'])):
                    for read in isoform['reads']:
                        cons[j] += read['data'][j]=='1'
                cons = [x/len(isoform['reads'])>0.3 for x in cons]
                if not True in cons:
                    continue
                record=list()
                record.append(tint['chrom'])
                record.append('freddie')
                record.append('transcript')
                f_seg = cons.index(True)
                record.append(str(tint['segs'][f_seg][0]))
                l_seg = len(cons)-1 - cons[::-1].index(True)
                record.append(str(tint['segs'][l_seg][1]))
                record.append('.')
                record.append('+')
                record.append('.')
                record.append('transcript_id "{chr}_{tint}_{iid}"; read_support "{read_support}";'.format(
                    chr=tint['chrom'],
                    tint=tint['id'],
                    iid=isoform['id'],
                    read_support=len(isoform['reads']),
                ))
                out_file.write('\t'.join(record))
                out_file.write('\n')
                eid = 1
                for d,group in groupby(enumerate(cons), lambda x: x[1]):
                    if d != True:
                        continue
                    group = list(group)
                    f_seg_idx = group[0][0]
                    l_seg_idx = group[-1][0]
                    record=list()
                    record.append(tint['chrom'])
                    record.append('freddie')
                    record.append('exon')
                    record.append(str(tint['segs'][f_seg][0]))
                    record.append(str(tint['segs'][l_seg][1]))
                    record.append('.')
                    record.append('transcript_id "{chr}_{tint}_{iid}"; exon_number "{eid}"; exon_id "{chr}_{tint}_{iid}_{eid}"; '.format(
                        chr=tint['chrom'],
                        tint=tint['id'],
                        iid=isoform['id'],
                        eid=eid,
                    ))
                    out_file.write('\t'.join(record))
                    out_file.write('\n')
                    eid+=1
    out_file.close()

def main():
    args = parse_args()

    tints = get_tints(args.cluster_tsv, args.segment_tsv)
    output_gtf(tints, args.output)



if __name__ == "__main__":
    main()
