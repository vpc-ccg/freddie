#!/usr/bin/env python3
import os
import argparse
import re
from itertools import groupby
from math import exp

smooth = list()
threshold = 0.9
while True:
    x = len(smooth)
    y = threshold/(1 + ((threshold-.5)/.5)*exp(-0.05*x))
    if x>5 and x*(threshold-y)<0.5:
        break
    smooth.append(round(y,2))
    assert len(smooth)<1000

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
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default='freddie_summarize',
                        help="Output file. Default: freddie_summarize.tsv")
    args = parser.parse_args()
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
        transcripts[tid]['tints'] = dict()
    return transcripts

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
            iid=line[5]
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
                iid=iid,
                data=data,
                poly_tail=line[6+1+len(data):],
            ))
    return tints

def process_tint(tint, transcripts):
    tint['tids'] = list()

    tids=set()
    tint_id =tint['id']
    for isoform in tint['isoforms'].values():
        for read in isoform['reads']:
            tid = read['tid']
            if not tid in transcripts:
                continue
            tids.add(tid)
            transcript=transcripts[tid]
            if not tint_id in transcript['tints']:
                 transcript['tints'][tint_id]=dict(
                    cov=0,
                    data=[0 for _ in tint['segs']],
                    missed=0,
                 )
            transcript['tints'][tint_id]['cov'] += 1
    pos_to_seg = list()
    for idx,(s,e) in enumerate(tint['segs']):
        for i in range(s,e):
            pos_to_seg.append(idx)
    tint_s,tint_e = tint['segs'][0][0],tint['segs'][-1][-1]
    for tid in tids:
        transcript=transcripts[tid]
        # if transcript['tints'][tint_id]['cov'] < 0:
        #     continue
        tint['tids'].append(tid)
        for s,e in transcript['intervals']:
            for i in range(s,e):
                if tint_s<=i<tint_e:
                    seg_idx = pos_to_seg[i-tint_s]
                    transcript['tints'][tint_id]['data'][seg_idx]+=1
                else:
                    transcript['tints'][tint_id]['missed']+=1
        for idx,(s,e) in enumerate(tint['segs']):
            l=e-s+1
            c=transcript['tints'][tint_id]['data'][idx]
            t=threshold
            if l < len(smooth):
                t = smooth[l]
            d = '2'
            if c/l > t:
                d = '1'
            elif c/l < 1-t:
                d = '0'
            transcript['tints'][tint_id]['data'][idx]=d
        transcript['tints'][tint_id]['data'] = ''.join(transcript['tints'][tint_id]['data'])

def process_overlaps(tint, transcripts):
    pos_to_seg = list()
    # for tid in tint['tids']:
    #     data = [0 for _ in tint['segs']]
    #     for s,e in transcripts[tid]['intervals']



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

def output_tint(tint, transcripts, outfile):
    record = ['tint']
    record.append(str(tint['id']))
    record.append(tint['chrom'])
    rid_count = sum(len(isoform['reads']) for isoform in tint['isoforms'].values())
    record.append(str(rid_count))
    tid_count = len(tint['tids'])
    record.append(str(tid_count))
    segs = ','.join(['{}-{}'.format(s,e) for s,e in tint['segs']])
    record.append(segs)
    print('\t'.join(record), file=outfile)

    for tid in tint['tids']:
        record = ['tid']
        transcript = transcripts[tid]
        record.append(str(tint['id']))
        record.append(transcript['chrom'])
        record.append(tid)
        record.append(transcript['name'])
        record.append(transcript['tints'][tint['id']]['data'])
        record.append(str(transcript['tints'][tint['id']]['cov']))
        record.append(str(transcript['tints'][tint['id']]['missed']))
        record.append(','.join(list(map(str,transcript['tints'].keys()))))
        print('\t'.join(record), file=outfile)

    for iid,isoform in tint['isoforms'].items():
        record = ['isoform']
        record.append(str(tint['id']))
        record.append(tint['chrom'])
        iname = '{}.{}'.format(tint['id'],iid)
        record.append(iname)
        record.append(isoform['data'])
        record.append(str(len(isoform['reads'])))
        print('\t'.join(record), file=outfile)
        for read in isoform['reads']:
            record = ['read']
            record.append(str(tint['id']))
            record.append(tint['chrom'])
            record.append(iname)
            record.append(read['name'])
            record.append(read['data'])
            print('\t'.join(record), file=outfile)

def main():
    args = parse_args()

    transcripts = get_transcripts(gtf=args.annotation_gtf)
    tints =  get_tints(cluster_tsv=args.cluster_tsv, segment_tsv=args.segment_tsv)
    for idx,tint in enumerate(tints.values()):
        process_tint(tint, transcripts)
    outfile = open(args.output, 'w+')
    for tint_id,tint in sorted(tints.items()):
        output_tint(tint, transcripts, outfile)
if __name__ == "__main__":
    main()
