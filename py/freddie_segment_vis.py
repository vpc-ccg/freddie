#!/usr/bin/env python3
import argparse
import re
import pickle
import numpy as np
from itertools import groupby

def parse_args():
    parser = argparse.ArgumentParser(
        description="Outputs pickle file for vis purposes")
    parser.add_argument("-s",
                        "--split-tsv",
                        type=str,
                        required=True,
                        help="Freddie split TSV file")
    parser.add_argument("-g",
                        "--segment-tsv",
                        type=str,
                        required=True,
                        help="Freddie segment TSV file")
    parser.add_argument("-a",
                        "--annotation-gtf",
                        type=str,
                        required=True,
                        help="Annotation GTF file path")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default='vis_segmentation.pickle',
                        help="Output path. Default: vis_segmentation.pickle")
    args = parser.parse_args()

    return args

# def read_annotation_segmentation(annotation_gtf, reads):
def read_annotation_gtf(annotation_gtf):
    cid_to_transcripts = dict()
    for line in open(annotation_gtf):
        if line[0]=='#':
            continue
        line = line.split('\t')
        if not line [2] == 'exon':
            continue
        chrom = line[0]
        if chrom not in cid_to_transcripts:
            cid_to_transcripts[chrom]=dict()
        gid = re.search(r'gene_id \"(?P<gid>ENSG\d{11})\"', line[8]).group('gid')
        tid = re.search(r'transcript_id \"(?P<tid>ENST\d{11})\"', line[8]).group('tid')
        if tid not in cid_to_transcripts[chrom]:
            cid_to_transcripts[chrom][tid]=dict(
                tid=tid,
                gid=gid,
                intervals=list(),
            )
        s,e=int(line[3]), int(line[4])
        cid_to_transcripts[chrom][tid]['intervals'].append((s,e))
    return cid_to_transcripts

def get_annotation_positions(cid_to_transcripts, w=5):
    cid_to_positions = dict()
    for cid,transcripts in cid_to_transcripts.items():
        pos_cnt = dict()
        for transcript in transcripts.values():
            for s,e in transcript['intervals']:
                pos_cnt[s] = pos_cnt.get(s,0)+1
                pos_cnt[e] = pos_cnt.get(e,0)+1
        positions = sorted(pos_cnt.keys())
        neighborhoods = [list()]
        idx = 0
        while idx < len(positions)-1:
            p_curr = positions[idx]
            p_next = positions[idx+1]
            if p_next-p_curr < w:
                if len(neighborhoods[-1]) == 0:
                    neighborhoods[-1].append(p_curr)
                else:
                    neighborhoods[-1][-1]==p_curr
                neighborhoods[-1].append(p_next)
            elif len(neighborhoods[-1]) != 0:
                neighborhoods.append(list())
            idx+=1
        if len(neighborhoods[-1])==0:
            neighborhoods.pop()
        final_positions = set(pos_cnt.keys())
        for neighborhood in neighborhoods:
            final_positions.difference_update(neighborhood)
            final_positions.add(int(round(np.average(neighborhood, weights=[pos_cnt[p] for p in neighborhood]))))
        final_positions = sorted(final_positions)
        for idx,p in enumerate(final_positions[:-1]):
            assert final_positions[idx+1]-p >= w
        cid_to_positions[cid]=final_positions
    return cid_to_positions

def get_segmentation_position(segment_tsv):
    cid_to_positions = dict()
    for line in open(segment_tsv):
        if line[0]!='#':
            continue
        line = line[1:].rstrip().split('\t')
        chrom = line[0]
        if not chrom in cid_to_positions:
            cid_to_positions[chrom]=set()
        cid_to_positions[chrom].update((int(x) for x in line[2].split(',')))
    return cid_to_positions

def switch_to_nearest(cid_to_s_pos, cid_to_a_pos, w=5):
    for chrom in cid_to_s_pos:
        cid_to_s_pos[chrom] = sorted(cid_to_s_pos[chrom])
        cid_to_s_pos[chrom] = [a for a,b in zip(cid_to_s_pos[chrom][:-1],cid_to_s_pos[chrom][1:]) if  b-a>w] + [cid_to_s_pos[chrom][-1]]
        positions = {p for p in cid_to_a_pos[chrom]}
        for idx,p in enumerate(cid_to_s_pos[chrom]):
            hits = [(abs(x-p),x) for x in range(p-w,p+w+1) if x in positions]
            if len(hits)!=0:
                cid_to_s_pos[chrom][idx]=min(hits)[1]

def buffer_print(positions, track, marks, size=200):
    buffer = ['','','']
    for p,t,m in zip(positions,track,marks):
        if max(map(len,buffer)) >= size:
            print('\n'.join(buffer))
            print()
            buffer = ['','','']
        if buffer[0]=='':
            buffer[0]=str(p)
        buffer[1]+=m
        buffer[2]+=str(t)
    if max(map(len,buffer)) > 0:
        print('\n'.join(buffer))

def get_seg_track(cid_to_s_pos, cid_to_a_pos):
    cid_to_segs = dict()
    for chrom in cid_to_s_pos:
        cid_to_segs[chrom]=dict(
                segs=list(),
                track=list()
        )
        final_positions = {0:3}
        for p in cid_to_s_pos[chrom]+cid_to_a_pos[chrom]:
            final_positions[p] = 0
        for p in cid_to_s_pos[chrom]:
            final_positions[p] |= 1
        for p in cid_to_a_pos[chrom]:
            final_positions[p] |= 2

        final_positions=sorted(final_positions.items())
        for (p1,t1),(p2,t2) in zip(final_positions[:-1],final_positions[1:]):
            cid_to_segs[chrom]['segs'].append((p1,p2))
            if t1 == 2 and t2 ==2:
                cid_to_segs[chrom]['track'].append('-')
            if t1 != 2 and t2 ==2:
                cid_to_segs[chrom]['track'].append('<')
            if t1 == 2 and t2 !=2:
                cid_to_segs[chrom]['track'].append('>')
        # # a_only_intervals = list()
        # # for d,group in groupby(enumerate(track), lambda x: x[1]):
        # #     if d != 2:
        # #         continue
        # #     group = list(group)
        # #     f_seg_idx = group[0][0]
        # #     l_seg_idx = group[-1][0]
        # #     a_only_intervals.append((f_seg_idx,l_seg_idx))
        # marks=[' ' if t !=2 else '-' for t in track]
        # # for f,l in a_only_intervals:
        # #     if f > 0:
        # #         marks[f-1]='<'
        # #     if l < len(track)-1:
        # #         marks[l+1]='>'
        # cid_to_positions[chrom]=dict(
        #     positions=positions,
        #     marks=marks,
        # )
    return cid_to_segs

def get_reads(split_tsv):
    cid_to_reads = dict()
    for line in open(split_tsv):
        if line[0]=='#':
            continue
        line=line.rstrip().split('\t')
        read=dict(
            rid=int(line[0]),
            name=line[1],
            tid=line[1].split('_')[0],
            strand=line[3],
            tint=line[4],
            intervals = list(),
        )
        chrom = line[2]
        for interval in line[5:]:
            interval = interval.split(':')[0].split('-')
            s = int(interval[0])
            e = int(interval[1])
            read['intervals'].append((s,e))
        if not chrom in cid_to_reads:
            cid_to_reads[chrom]=list()
        cid_to_reads[chrom].append(read)
    return cid_to_reads

def get_data(intervals, segs):
    locs = set()
    for s,e in intervals:
        locs.update(range(s,e))
    f_loc = min(locs)
    l_loc = max(locs)

    data = dict()
    for idx,(s,e) in enumerate(segs):
        flag = False
        for s2,e2 in intervals:
            if s<=s2<=e or s2<=s<=e2:
                flag=True
                break
        if flag == False:
            continue
        c=sum(1 if p in locs else 0 for p in range(s,e))/(e-s)
        if c > 0.9:
            data[idx]=1
        elif c < 0.1:
            data[idx]=0
        else:
            data[idx]=2
    return data

def main():
    args = parse_args()
    cid_to_transcripts = read_annotation_gtf(args.annotation_gtf)
    cid_to_a_pos = get_annotation_positions(cid_to_transcripts)
    cid_to_s_pos = get_segmentation_position(args.segment_tsv)
    switch_to_nearest(cid_to_s_pos, cid_to_s_pos)
    cid_to_segs=get_seg_track(cid_to_s_pos, cid_to_s_pos)
    cid_to_reads = get_reads(args.split_tsv)

    for chrom,reads in cid_to_reads.items():
        for idx,read in enumerate(reads):
            if idx%500==0:
                print('Chrom {}: Read {}/{}'.format(chrom,idx,len(reads)))
            read['data']=get_data(intervals=read['intervals'], segs=cid_to_segs[chrom]['segs'])
        for idx,transcript in enumerate(cid_to_transcripts[chrom].values()):
            if idx%500==0:
                print('Chrom {}: Transcript {}/{}'.format(chrom,idx,len(cid_to_transcripts[chrom])))
            transcript['data']=get_data(intervals=transcript['intervals'], segs=cid_to_segs[chrom]['segs'])

    with open(args.output, 'wb+') as out_file:
        pickle.dump((cid_to_segs,cid_to_transcripts,cid_to_reads), out_file)

if __name__ == "__main__":
    main()
