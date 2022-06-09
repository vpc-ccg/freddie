#!/usr/bin/env python3
import os
import glob
from multiprocessing import Pool

import argparse
from itertools import groupby


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract alignment information from BAM/SAM file and splits reads into distinct transcriptional intervals")
    parser.add_argument("-s",
                        "--split-dir",
                        type=str,
                        required=True,
                        help="Path to directory of Freddie segment")
    parser.add_argument("-c",
                        "--cluster-dir",
                        type=str,
                        required=True,
                        help="Path to directory of Freddie cluster")
    parser.add_argument("-m",
                        "--majority-threshold",
                        type=float,
                        default=0.50,
                        help="Majority threshold of reads to adjust exon boundary using the original alignments. Default: 0.5")
    parser.add_argument("-w",
                        "--correction-window",
                        type=int,
                        default=8,
                        help="The +/- window around segment boundary to look for read alignment boundaries to use for correcting the exon boundaries. Value of 0 means no correction. Default 8.")
    parser.add_argument("-t",
                        "--threads",
                        type=int,
                        default=1,
                        help="Number of threads to use")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default='freddie_isoforms.gtf',
                        help="Path to output file. Default: freddie_isoforms.gtf")
    args = parser.parse_args()
    assert 0.5 <= args.majority_threshold <= 1.0
    assert 0 <= args.correction_window <= 20
    assert 0 < args.threads
    return args


def run_consensus(consensus_args):
    (
        contig,
        tint_id,
        cluster_tsv,
        split_tsv,
        majority_threshold,
        correction_window,
    ) = consensus_args
    print(f'Building isoforms for contig {contig}')
    segments, reads, isoforms = read_cluster(cluster_tsv)
    isoforms_cons(isoforms, segments, reads)
    read_split(split_tsv, reads)

    correct_boundaries('starts', isoforms, reads,
                       majority_threshold, correction_window)
    correct_boundaries('ends', isoforms, reads,
                       majority_threshold, correction_window)

    return get_gtf_records(isoforms)


def get_gtf_records(isoforms):
    gtf_records = list()
    for isoform_key, isoform in isoforms.items():
        isoform_record = list()
        if not 'starts' in isoform:
            continue
        chrom, tint, pid, iid = isoform_key
        starts = isoform['starts']
        ends = isoform['ends']
        strand = isoform['strand']
        transcript_name = '{chrom}_{tint}_{iid}'.format(
            chrom=chrom,
            tint=tint,
            iid=iid,
        )
        gtf_record_key = (chrom, starts[0])

        record = list()
        record.append(chrom)
        record.append('freddie')
        record.append('transcript')
        record.append(str(starts[0]+1))
        record.append(str(ends[-1]))
        record.append('.')
        record.append(strand)
        record.append('.')
        record.append('transcript_id "{transcript_name}"; read_support "{read_support}";'.format(
            transcript_name=transcript_name,
            read_support=len(isoform['rids']),
        ))
        isoform_record.append('\t'.join(record))
        for eid, (s, e) in enumerate(zip(starts, ends), start=1):
            record = list()
            record.append(chrom)
            record.append('freddie')
            record.append('exon')
            record.append(str(s))
            record.append(str(e))
            record.append('.')
            record.append(strand)
            record.append('.')
            record.append('transcript_id "{transcript_name}"; exon_number "{eid}"; exon_id "{transcript_name}_{eid}"; '.format(
                transcript_name=transcript_name,
                eid=eid,
            ))
            isoform_record.append('\t'.join(record))
        gtf_records.append((gtf_record_key, '\n'.join(isoform_record)))
    return gtf_records


def correct_boundaries(side, isoforms, reads, majority_threshold, correction_window):
    if correction_window == 0:
        return
    assert side in ['starts', 'ends']
    for isoform in isoforms.values():
        if not side in isoform:
            continue
        for idx, iso_s in enumerate(isoform[side]):
            cur = {
                x: 0 for x in range(-correction_window, correction_window+1)}
            for rid in isoform['rids']:
                for read_s in reads[rid][side]:
                    x = read_s - iso_s
                    if not x in cur:
                        continue
                    cur[x] += 1
            for x, v in cur.items():
                if v/len(isoform['rids']) >= majority_threshold:
                    isoform[side][idx] = x+iso_s


def read_split(split_tsv, reads):
    for line in open(split_tsv):
        if line.startswith('#'):
            continue
        line = line.rstrip().split('\t')
        rid = int(line[0])
        if not rid in reads:
            continue
        intervals = [i.split(':')[0].split('-') for i in line[5:]]
        starts, ends = zip(*[(int(i[0]), int(i[1])) for i in intervals])
        reads[rid]['starts'] = starts
        reads[rid]['ends'] = ends
        for s, e in zip(starts, ends):
            assert s < e


def read_cluster(cluster_tsv):
    segments = dict()
    reads = dict()
    isoforms = dict()
    for line in open(cluster_tsv):
        line = line.rstrip().split('\t')
        if line[0][0] == '#':
            chrom = line[0][1:]
            tint = int(line[1])
            segments[(chrom, tint)] = [int(x) for x in line[2].split(',')]
            segments[(chrom, tint)] = [(s, e) for s, e in zip(
                segments[(chrom, tint)][:-1], segments[(chrom, tint)][1:])]
            continue
        if line[0].startswith('isoform_'):
            continue
        if line[7] == '*':
            continue
        read = dict()
        read['rid'] = int(line[0])
        read['rname'] = line[1]
        read['chrom'] = line[2]
        read['strand'] = line[3]
        read['tint'] = int(line[4])
        read['pid'] = int(line[5])
        read['tail'] = line[6]
        read['iid'] = int(line[7])
        read['data'] = line[8]
        assert len(read['data']) == len(
            segments[(read['chrom'], read['tint'])])
        reads[read['rid']] = read
        isoform_key = (read['chrom'], read['tint'], read['pid'], read['iid'])
        if not isoform_key in isoforms:
            isoforms[isoform_key] = dict(rids=set())
        isoforms[isoform_key]['rids'].add(read['rid'])

    for isoform_key, isoform in isoforms.items():
        l = set()
        for rid in isoform['rids']:
            l.add(len(reads[rid]['data']))
        assert len(l) == 1

    return segments, reads, isoforms


def isoforms_cons(isoforms, segments, reads):
    for isoform_key, isoform in isoforms.items():
        chrom, tint, _, _ = isoform_key
        cons = [0 for _ in segments[(chrom, tint)]]
        cov = [0 for _ in segments[(chrom, tint)]]
        M = len(segments[(chrom, tint)])
        N = len(isoform['rids'])
        tails = {'N': 0, 'S': 0, 'E': 0}
        for rid in isoform['rids']:
            read = reads[rid]
            assert len(read['data']) == M, (M, isoform_key, read)
            if not '1' in read['data']:
                continue
            if read['tail'] == 'S':
                first = 0
            else:
                first = read['data'].index('1')
            if read['tail'] == 'S':
                last = M - 1
            else:
                last = M - 1 - read['data'][::-1].index('1')
            assert 0<=first<=last<M,f'first {first}, last {last}, M {M}'
            for j in range(first,last+1):
                X = read['data'][j] == '1'
                cons[j] += X
                cov[j] += 1
            tails[read['tail']] += 1
        
        cons = [x/c > 0.5 if x >= 3 else False for x,c in zip(cons,cov) ]
        if not True in cons:
            continue
        if tails['S'] > tails['E']:
            isoform['strand'] = '-'
        else:
            isoform['strand'] = '+'
        starts = list()
        ends = list()
        for d, group in groupby(enumerate(cons), lambda x: x[1]):
            if d != True:
                continue
            group = list(group)
            f_seg_idx = group[0][0]
            l_seg_idx = group[-1][0]
            starts.append(segments[(chrom, tint)][f_seg_idx][0])
            ends.append(segments[(chrom, tint)][l_seg_idx][1])
        isoform['starts'], isoform['ends'] = starts, ends
        for s, e in zip(starts, ends):
            assert s < e, (s, e)


def main():
    args = parse_args()

    consensus_args = list()
    for contig in os.listdir(args.cluster_dir):
        if not os.path.isdir('{}/{}'.format(args.cluster_dir, contig)):
            continue
        for cluster_tsv in glob.iglob('{}/{}/cluster_*.tsv'.format(args.cluster_dir, contig)):
            tint_id = int(cluster_tsv[:-4].split('/')[-1].split('_')[-1])
            split_tsv = '{}/{}/split_{}_{}.tsv'.format(args.split_dir, contig, contig, tint_id)
            assert os.path.isfile(split_tsv), split_tsv
            consensus_args.append([
                contig,
                tint_id,
                cluster_tsv,
                split_tsv,
                args.majority_threshold,
                args.correction_window,
            ])

    gtf_records = list()
    if args.threads > 1:
        p = Pool(args.threads)
        for idx, tint_gtf_records in enumerate(p.imap_unordered(run_consensus, consensus_args, chunksize=5)):
            gtf_records.extend(tint_gtf_records)
    else:
        for idx, tint_gtf_records in enumerate(map(run_consensus, consensus_args)):
            gtf_records.extend(tint_gtf_records)
    gtf_records.sort()

    outfile = open(args.output, 'w+')
    for (k,record) in gtf_records:
        outfile.write(record)
        outfile.write('\n')
    outfile.close()


if __name__ == "__main__":
    main()
