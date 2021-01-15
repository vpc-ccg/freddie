#!/usr/bin/env python3
import argparse
import glob
import random

def parse_args():
    parser = argparse.ArgumentParser(
        description="Converts GTF to BED files")
    parser.add_argument("-s",
                        "--seqpare",
                        type=str,
                        required=True,
                        help="Path to seqpare scores TSV file")
    parser.add_argument("-b",
                        "--beds",
                        type=str,
                        required=True,
                        help="Path to BED directory of the tool")
    parser.add_argument("-r",
                        "--ref-beds",
                        type=str,
                        required=True,
                        help="Path to BED directory of the reference")
    parser.add_argument("-t",
                        "--trials",
                        type=int,
                        default=100,
                        help="Number of trials. Default: 100")
    parser.add_argument("-m",
                        "--min-seqpare",
                        type=float,
                        default=0.6,
                        help="Min acceptable seqpare score. Default: 0.5")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default='freddie_pair.tsv',
                        help="Path to output. Default: freddie_pair.tsv")
    args = parser.parse_args()
    assert args.trials > 0
    assert 1.0 >= args.min_seqpare >= 0.0
    args.beds = args.beds.rstrip('/')
    args.ref_beds = args.ref_beds.rstrip('/')
    return args

def main():
    args = parse_args()

    random.seed(42)
    iids = [f.split('/')[-1][:-4] for f in glob.glob('{}/*.bed'.format(args.beds))]
    tids = [f.split('/')[-1][:-4] for f in glob.glob('{}/*.bed'.format(args.ref_beds))]
    seqpare =  [l.split('\t') for l in open(args.seqpare) if float(l.split('\t')[2]) >= args.min_seqpare]
    pairings = dict()
    for idx in range(args.trials):
        iid_to_tid = {iid:'None' for iid in iids}
        tid_to_iid = {tid:'None' for tid in tids}
        scores = sorted([(float(x[2]),random.random(),x[0],x[1]) for x in seqpare], reverse=True)
        for s,_,tid,iid in scores:
            if iid_to_tid[iid] != 'None':
                continue
            if tid_to_iid[tid] != 'None':
                continue
            iid_to_tid[iid]=tid
            tid_to_iid[tid]=iid
            if not (iid,tid) in pairings:
                pairings[(iid,tid)] = dict(
                    s=s,
                    c=list(),
                )
            assert pairings[(iid,tid)]['s']==s
            pairings[(iid,tid)]['c'].append(idx)
        for iid,tid in iid_to_tid.items():
            if not tid == 'None':
                continue
            if not (iid,tid) in pairings:
                pairings[(iid,tid)] = dict(
                    s=0.0,
                    c=list(),
                )
            pairings[(iid,tid)]['c'].append(idx)
        for tid,iid in tid_to_iid.items():
            if not iid == 'None':
                continue
            if not (iid,tid) in pairings:
                pairings[(iid,tid)] = dict(
                    s=0.0,
                    c=list(),
                )
            pairings[(iid,tid)]['c'].append(idx)

    pairings = sorted([(v['c'], v['s'], k[0],k[1]) for k,v in pairings.items()], reverse=True)
    out_file = open(args.output, 'w+')
    for c,s,iid,tid in pairings:
        print(c,s,iid,tid, sep='\t', file=out_file)
    out_file.close()

if __name__ == "__main__":
    main()
