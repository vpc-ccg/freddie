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
    parser.add_argument("-p",
                        "--split-tsv",
                        type=str,
                        required=True,
                        help="Path to TSV file of Freddie split")
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
                        default='freddie_summarize.tsv',
                        help="Output file. Default: freddie_summarize.tsv")
    parser.add_argument("-ob",
                        "--beds-dir",
                        type=str,
                        default='freddie_beds/',
                        help="Output directory. Will be created if does not exist. Default: freddie_beds/")
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

def get_tints(split_tsv, cluster_tsv, segment_tsv):
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
            assert tint_id in tints
            tints[tint_id]['isoforms'][iid]=dict(
                id=iid,
                data=data,
                reads=list(),
            )
        else:
            line=line.rstrip().split('\t')
            rid=int(line[0])
            tint_id=int(line[4])
            iid=line[7]
            if iid =='*':
                iid = 'garbage'
            data=rid_to_data[rid]
            assert tint_id in tints
            assert iid in tints[tint_id]['isoforms'], (line,tint_id,iid,tints[tint_id]['isoforms'].keys())
            tints[tint_id]['isoforms'][iid]['reads'].append(dict(
                rid=rid,
                name=line[1],
                tid=line[1].split('_')[0],
                chrom=line[2],
                strand=line[3],
                tint=tint_id,
                partition=int(line[5]),
                poly_tail_category=line[6],
                iid=iid,
                data=data,
                poly_tail=line[8+1+len(data):],
            ))
    cinterval_re   = '([0-9]+)-([0-9]+)'
    tint_prog   = re.compile(r'#%(chr_re)s\t%(cid_re)s\t%(intervals_re)s\t%(read_count_re)s\n$' % {
        'chr_re'        : '(?P<chr>[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)',
        'cid_re'        : '(?P<cid>[0-9]+)',
        'intervals_re'  : '(?P<intervals>%(i)s(,%(i)s)*)' % {'i':cinterval_re},
        'read_count_re' : '(?P<read_count>[0-9]+)',
    })
    cinterval_prog  = re.compile(cinterval_re)
    for line in open(split_tsv):
        if not line[0]=='#':
            continue
        re_dict = tint_prog.match(line).groupdict()
        tint_id = int(re_dict['cid'])
        tints[tint_id]['intervals']  = [(int(x[0]),int(x[1])) for x in cinterval_prog.findall(re_dict['intervals'])]
    return tints

def process_tint(tint, transcripts):
    tint['tids'] = list()
    tids=dict()

    tint_id =tint['id']
    for isoform in tint['isoforms'].values():
        for read in isoform['reads']:
            tid = read['tid']
            if not tid in transcripts:
                continue
            transcript=transcripts[tid]
            if not tid in tids:
                tids[tid]=dict(
                    read_cov=0,
                    data=[0 for _ in tint['segs']],
                    covered=0,
                    missed=0,
                 )
            tids[tid]['read_cov'] += 1
    pos_to_seg = dict()
    tint_skips = {(e1,s2):False for (s1,e1),(s2,e2) in zip(tint['intervals'][:-1],tint['intervals'][1:])}
    for idx,(s,e) in enumerate(tint['segs']):
        if (s,e) in tint_skips:
            assert tint_skips[(s,e)]==False
            tint_skips[(s,e)]=True
        for i in range(s,e):
            pos_to_seg[i]=idx
    assert all(tint_skips.values()),(tint_skips,tint['intervals'])

    for tid in tids:
        transcript=transcripts[tid]
        if transcript['chrom'] != tint['chrom']:
            continue
        for s,e in transcript['intervals']:
            for i in range(s,e):
                if i in pos_to_seg:
                    seg_idx = pos_to_seg[i]
                    tids[tid]['data'][seg_idx]+=1
                    tids[tid]['covered']+=1
                else:
                    tids[tid]['missed']+=1
        for idx,(s,e) in enumerate(tint['segs']):
            l=e-s+1
            c=tids[tid]['data'][idx]
            t=threshold
            if l < len(smooth):
                t = smooth[l]
            d = '2'
            if c/l > t:
                d = '1'
            elif c/l < 1-t:
                d = '0'
            tids[tid]['data'][idx]=d
        tids[tid]['data'] = ''.join(tids[tid]['data'])
    for tid in tids:
        if tids[tid]['covered'] < 20:
            continue
        tint['tids'].append(tid)
        transcripts[tid]['tints'][tint['id']]=tids[tid]

def output_beds(tint, transcripts, tid_dir, iso_dir):
    os.makedirs(tid_dir, exist_ok=True)
    os.makedirs(iso_dir, exist_ok=True)

    chrom = tint['chrom']
    for tid in tint['tids']:
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
        record.append(str(transcript['tints'][tint['id']]['read_cov']))
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
    tints =  get_tints(split_tsv=args.split_tsv, cluster_tsv=args.cluster_tsv, segment_tsv=args.segment_tsv)
    for idx,tint in enumerate(tints.values()):
        process_tint(tint, transcripts)
    outfile = open(args.output, 'w+')
    for tint_id,tint in sorted(tints.items()):
        output_tint(
            tint        = tint,
            transcripts = transcripts,
            outfile     = outfile,
        )
        output_beds(
            tint        = tint,
            transcripts = transcripts,
            tid_dir     = '{}/{}/tids'.format(args.beds_dir, tint['id']),
            iso_dir     = '{}/{}/isos'.format(args.beds_dir, tint['id']),
        )
if __name__ == "__main__":
    main()
