#!/usr/bin/env python3
import argparse
from itertools import groupby

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract alignment information from BAM/SAM file and splits reads into distinct transcriptional intervals")
    parser.add_argument("-c",
                        "--cluster-tsv",
                        type=str,
                        required=True,
                        help="Path to Freddie cluster TSV file of the reads")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default='freddie_split.tsv',
                        help="Path to output file. Default: freddie_isoforms.gtf")
    args = parser.parse_args()
    return args

def read_cluster(cluster_tsv):
    tints = dict()
    for line in open(cluster_tsv):
        if line.startswith('#'):
            chr,tint_id,segs=line[1:-1].split('\t')
            tint_id = int(tint_id)
            segs = segs.split(',')
            segs = [(int(s),int(e)) for s,e in zip(segs[:-1], segs[1:])]
            tints[tint_id]=dict(
                id=tint_id,
                chr=chr,
                segs=segs,
                isoforms=dict()
            )
        elif line.startswith('isoform_'):
            iid,tint_id,data = line[8:-1].split('\t')
            tint_id=int(tint_id)
            iid=int(iid)
            tints[tint_id]['isoforms'][iid]=dict(
                id=iid,
                data=data,
            )
    return tints

def output_gtf(tints, outpath):
    out_file=open(outpath, 'w+')
    for tint in tints.values():
        for isoform in tint['isoforms'].values():
            if not '1' in isoform['data']:
                continue
            record=list()
            record.append(tint['chr'])
            record.append('freddie')
            record.append('transcript')
            f_seg = isoform['data'].index('1')
            record.append(str(tint['segs'][f_seg][0]))
            l_seg = len(isoform['data'])-1 - isoform['data'][::-1].index('1')
            record.append(str(tint['segs'][l_seg][1]))
            record.append('.')
            record.append('+')
            record.append('.')
            record.append('transcript_id "{chr}_{tint}_{iid}"; '.format(
                chr=tint['chr'],
                tint=tint['id'],
                iid=isoform['id'],
            ))
            out_file.write('\t'.join(record))
            out_file.write('\n')
            eid = 1
            for d,group in groupby(enumerate(isoform['data']), lambda x: x[1]):
                if d != '1':
                    continue
                group = list(group)
                f_seg_idx = group[0][0]
                l_seg_idx = group[-1][0]
                record=list()
                record.append(tint['chr'])
                record.append('freddie')
                record.append('exon')
                record.append(str(tint['segs'][f_seg][0]))
                record.append(str(tint['segs'][l_seg][1]))
                record.append('.')
                record.append('transcript_id "{chr}_{tint}_{iid}"; exon_number "{eid}"; exon_id "{chr}_{tint}_{iid}_{eid}"; '.format(
                    chr=tint['chr'],
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

    tints = read_cluster(args.cluster_tsv)
    output_gtf(tints, args.output)



if __name__ == "__main__":
    main()
