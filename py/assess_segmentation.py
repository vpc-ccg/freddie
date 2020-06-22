#!/usr/bin/env python3
import argparse
import subprocess
import re

def parse_args():
    parser = argparse.ArgumentParser(
        description="Outputs SAM alignments using deSALT splice aligner")
    # parser.add_argument("-r",
    #                     "--reads",
    #                     nargs="+",
    #                     type=str,
    #                     required=True,
    #                     help="Space separated paths to reads in FASTQ or FASTA format")
    parser.add_argument("-a",
                        "--annotation-gtf",
                        type=str,
                        required=True,
                        help="Annotation GTF file path")
    parser.add_argument("-s",
                        "--segment-tsv",
                        type=str,
                        required=True,
                        help="Freddie segmentation TSV file path")
    # parser.add_argument("-i",
    #                     "--isoform-gtf",
    #                     type=str,
    #                     required=True,
    #                     help="Tool isoform GTF path")
    parser.add_argument("-t",
                        "--tid-cov-threshold",
                        type=int,
                        default=4,
                        help="Min read coverage threshold per transcript. Default: 4")
    args = parser.parse_args()

    return args

# def read_annotation_segmentation(annotation_gtf, reads):
def read_annotation_gtf(annotation_gtf, tid_cov, threshold=4):
    positions = dict()
    f_pos = dict()
    l_pos = dict()
    for line in open(annotation_gtf):
        if line[0]=='#':
            continue
        line = line.split('\t')
        if not line [2] == 'exon':
            continue
        tid = re.search(r'transcript_id \"(?P<tid>ENST\d{11})\"', line[8]).group('tid')
        if tid_cov.get(tid, 0) < threshold:
            continue
        chrom = line[0]
        if not chrom in positions:
            positions[chrom] = set()
        f_pos[tid] = min(int(line[3]), f_pos.get(tid, float('inf')))
        l_pos[tid] = max(int(line[4]), f_pos.get(tid, 0))
        for p in [int(line[3]), int(line[4])]:
            positions[chrom].add(p)
    b_pos = set(f_pos.values()) | set(l_pos.values())
    return positions,b_pos

def read_isoform_segmentation(isoform_gtf):
    isoform_positions = dict()
    for line in open(isoform_gtf):
        line = line.split('\t')
        if not line[2] == 'exon':
            continue
        chrom = line[0]
        if not chrom in isoform_positions:
            isoform_positions[chrom] = set()
        isoform_positions[chrom].add(int(line[3]))
        isoform_positions[chrom].add(int(line[4]))
    # for k,v in isoform_positions.items():
    #     isoform_positions[k]=sorted(v)
    return isoform_positions

def read_segmentation_tsv(segment_tsv):
    positions = dict()
    tid_cov = dict()
    for line in open(segment_tsv):
        if line[0] == '#':
            line = line[1:].split('\t')
            chrom = line[0]
            if not chrom in positions:
                positions[chrom] = set()
            for p in line[2].split(','):
                positions[chrom].add(int(p))
        else:
            line = line.split('\t')
            tid=line[1].split('_')[0]
            tid_cov[tid]= tid_cov.get(tid, 0) + 1

    return positions,tid_cov

def compare(a_pos, b_pos, a_bound=set(), w=5):
    sums = dict()
    for chrom,positions in a_pos.items():
        for p in positions:
            s=sum(x in b_pos[chrom] for x in range(p-w,p+w+1))
            if p in a_bound:
                sums[(s,'boundary')] = sums.get((s,'boundary'), 0) + 1
            sums[(s,)] = sums.get((s,), 0) + 1

    return sums

    #     a_idx = 0
    #     b_idx = 0
    #     while a_idx < len(a_pos[chrom]) and b_idx < len(b_pos[chrom]):
    #         if a_pos[chrom][a_idx] < b_pos[chrom][b_idx]:
    #             left = a_pos[chrom][a_idx]
    #             right = ''
    #             a_idx += 1
    #         elif a_pos[chrom][a_idx] > b_pos[chrom][b_idx]:
    #             left = ''
    #             right = b_pos[chrom][b_idx]
    #             b_idx += 1
    #         else:
    #             left = a_pos[chrom][a_idx]
    #             right = b_pos[chrom][b_idx]
    #             a_idx += 1
    #             b_idx += 1
    #         print('{:>8} {}'.format(left,right))
    #     for v in a_pos[chrom][a_idx:]:
    #         print('{:>8} {}'.format(v,''))
    #     for v in b_pos[chrom][b_idx:]:
    #         print('{:>8} {}'.format('',v))
def main():
    args = parse_args()

    freddie_positions,tid_cov        = read_segmentation_tsv(segment_tsv=args.segment_tsv)
    annotation_positions,b_pos = read_annotation_gtf(annotation_gtf=args.annotation_gtf, tid_cov=tid_cov, threshold=args.tid_cov_threshold)
    Freddie_sum=sum(map(len, freddie_positions.values()))
    print('# of positions in Freddie:', Freddie_sum)

    annotation_sum=sum(map(len, annotation_positions.values()))
    print('# of positions in annotation: {} ({}% of which are boundary)'.format(annotation_sum, len(b_pos)*1000//annotation_sum/10.0))


    print('--> Freddie matches in annotation:')
    d=compare(a_pos=freddie_positions, b_pos=annotation_positions)
    for k,v in sorted(d.items()):
        print('{} ({}%) positions with {} matches'.format(v,v*1000//Freddie_sum/10.0,k))
    print('--> Annotation matches in freddie:')
    d=compare(a_pos=annotation_positions, b_pos=freddie_positions, a_bound=b_pos)
    for k,v in sorted(d.items()):
        print('{} ({}%) positions with {} matches'.format(v,v*1000//annotation_sum/10.0,k))



if __name__ == "__main__":
    main()
