#!/usr/bin/env python3
import os
import argparse
import re
from itertools import groupby

def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert cluster file into BED files and include ground truth overlapping GTF isoforms.")
    parser.add_argument("-a",
                        "--annotation-gtf",
                        type=str,
                        required=True,
                        help="Path to GTF file of annotations")
    parser.add_argument("-t",
                        "--tid-tint",
                        type=str,
                        required=True,
                        help="Path to TXT with TID and TINT ID in each line")
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
    parser.add_argument("-od",
                        "--out-dir",
                        type=str,
                        default='freddie_seqpare',
                        help="Output directory. Will be created if does not exist. Default: freddie_seqpare")
    args = parser.parse_args()
    args.out_dir = args.out_dir.rstrip('/')
    return args

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

def get_tints(cluster_tsv, segment_tsv, tid_tint_tsv):
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
                intersecting_tids=set(),
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
            iid=line[7]
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
                partition=int(line[5]),
                poly_tail_category=line[6],
                iid=iid,
                data=data,
                poly_tail=line[8+1+len(data):],
            ))
    for line in open(tid_tint_tsv):
        tid,tint = line.rstrip().split()
        tint = int(tint)
        tints[tint]['intersecting_tids'].add(tid)
    # for tint_id in tints.keys():
    #     tints[tint_id]['tids'] = set()
    #     for isoform in tints[tint_id]['isoforms'].values():
    #         for read in isoform['reads']:
    #             tints[tint_id]['tids'].add(read['tid'])
    return tints

def add_intersecting_tids(tint, transcripts):
    tint['intersecting_tids'] = list()

    tint_s,tint_e = tint['segs'][0][0],tint['segs'][-1][-1]
    for tid,transcript in transcripts.items():
        if transcript['chrom']!=tint['chrom']:
            continue
        s,e = transcript['intervals'][0][0],transcript['intervals'][-1][-1]
        if not (tint_e<s or e<tint_s):
            tint['intersecting_tids'].append(tid)
    tint['intersecting_tids'].sort(key=lambda tid: transcripts.get(tid, {'name':tid})['name'])

def output_beds(tint, transcripts, tid_dir, iso_dir):
    os.makedirs(tid_dir, exist_ok=True)
    os.makedirs(iso_dir, exist_ok=True)

    chrom = tint['chrom']
    for tid in tint['intersecting_tids']:
        name = transcripts[tid]['name']
        bed_file = open('{}/{}.bed'.format(tid_dir,name), 'w+')
        for s,e in transcripts[tid]['intervals']:
            print('{}\t{}\t{}'.format(
                chrom,
                s,
                e,
            ), file=bed_file)
        bed_file.close()
    for isoform in tint['isoforms'].values():
        if isoform['id']=='garbage':
            continue
        bed_file = open('{}/{}.bed'.format(iso_dir,isoform['id']), 'w+')
        cons = [sum(r['data'][i]=='1' for r in isoform['reads'])/len(isoform['reads'])>0.3 for i in range(len(tint['segs']))]
        for d,group in groupby(enumerate(cons), lambda x: x[1]):
            if d != True:
                continue
            group = list(group)
            f_seg_idx = group[0][0]
            l_seg_idx = group[-1][0]
            print('{}\t{}\t{}'.format(
                chrom,
                tint['segs'][f_seg_idx][0],
                tint['segs'][l_seg_idx][1],
            ), file=bed_file)
        bed_file.close()

def main():
    args = parse_args()

    transcripts = get_transcripts(gtf=args.annotation_gtf)
    # for tid,transcript in transcripts.items():
    #     print(tid)
    #     print(transcript['intervals'])
    #     print(transcript['name'])
    tints =  get_tints(cluster_tsv=args.cluster_tsv, segment_tsv=args.segment_tsv, tid_tint_tsv=args.tid_tint)
    for tint in tints.values():
        # add_intersecting_tids(tint, transcripts)
        output_beds(
            tint=tint,
            transcripts=transcripts,
            tid_dir='{}/{}/tids'.format(args.out_dir, tint['id']),
            iso_dir='{}/{}/isos'.format(args.out_dir, tint['id']),
        )
        # print(tint['id'])
        # print(tint['chrom'])
        # print(tint['segs'])
        # print(tint['tids'])
        # print(tint['intersecting_tids'])
        # print(tint.keys())
        # # for isoform in tint['isoforms'].values():
        # #     print(isoform['id'])
        # #     print(isoform['data'])
        # #     print(isoform.keys())
        # #     for read in isoform['reads']:
        # #         print(read)

if __name__ == "__main__":
    main()
