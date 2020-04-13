#!/usr/bin/env python3
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(
        description="Output PAF from SAM")
    parser.add_argument("-s",
                        "--sam",
                        type=str,
                        required=True,
                        help="Path to input SAM file")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to output PAF file")
    args = parser.parse_args()
    return args

def output_paf_from_sam(sam, paf):
    lengths = dict()
    outfile = open(paf, 'w+')
    for line in open(sam):
        line = line.rstrip().split('\t')
        if line[0][0]=='@':
            if line[0]=='@SQ':
                for f in line[1:]:
                    f = f.split(':')
                    if f[0]=='SN':
                        tname = f[1]
                    if f[0]=='LN':
                        tlen = int(f[1])
                lengths[tname] = tlen
            continue
        flag = int(line[1])
        if not flag < 256:
            continue
        if line[3] == '*':
            continue
        assert sum([len(x[0]+x[1]) for x in re.findall(r'(\d+)([M|I|D|N|S|H|P|=|X]{1})', line[5])])==len(line[5])
        cigar = [(int(x[0]),x[1]) for x in re.findall(r'(\d+)([M|I|D|N|S|H|P|=|X]{1})', line[5])]
        if len(cigar) == 0:
            continue
        qname = line[0]
        qstart = 0
        if cigar[0][1] in ['S']:
            qstart += cigar[0][0]
        qlen = 0
        for c,t in cigar:
            if t in ['I','S','M','=','X']:
                qlen+=c
        assert qlen==len(line[9])==len(line[10])
        qend = qlen
        if cigar[-1][1] in ['S']:
            qend -= cigar[-1][0]
        if flag & 16:
            strand = '-'
        else:
            strand = '+'
        tname = line[2]
        tlen = lengths[tname]
        tstart = int(line[3])-1
        tend = tstart + 1
        for c,t in cigar:
            if t in ['D','M','=','X']:
                tend+=c
        r_matches = 0

        qstart_c = qstart
        qend_c   = qstart
        tstart_c = tstart
        tend_c   = tstart

        intervals = list()
        for c,t in cigar:
            if t in ['D']:
                tend_c+=c
            if t in ['I']:
                qend_c+=c
            if t in ['X', '=', 'M']:
                tend_c+=c
                qend_c+=c
            if t == 'N':
                intervals.append((
                    qstart_c,
                    qend_c,
                    tstart_c,
                    tend_c,
                ))
                tend_c+=c
                tstart_c=tend_c
                qstart_c=qend_c
        if tstart_c < tend_c:
            intervals.append((
                qstart_c,
                qend_c,
                tstart_c,
                tend_c,
            ))
        for qstart_c,qend_c,tstart_c,tend_c in intervals:
            if strand == '-':
                qstart_c = qlen - qstart_c -1
                qend_c   = qlen - qend_c -1
                qstart_c,qend_c = qend_c,qstart_c
            print('\t'.join([
                str(qname),
                str(qlen),
                str(qstart_c),
                str(qend_c),
                str(strand),
                str(tname),
                str(tlen),
                str(tstart_c),
                str(tend_c),
                str(0),
                str(0),
                str(0),
                'tp:A:p',
                'oc:c:1'
            ]), file=outfile)
    outfile.close()

def main():
    args = parse_args()
    output_paf_from_sam(sam=args.sam, paf=args.paf)

if __name__ == "__main__":
    main()
