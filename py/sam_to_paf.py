#!/usr/bin/env python3
import argparse
import re
from itertools import zip_longest

cigar_re = re.compile(r'(\d+)([M|I|D|N|S|H|P|=|X]{1})')
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
        if line[2] == '*':
            continue
        assert sum([len(x[0]+x[1]) for x in cigar_re.findall(line[5])])==len(line[5]),'Something wrong with line:\n{}'.format('\n'.join(line))
        cigar = [(int(x[0]),x[1]) for x in cigar_re.findall(line[5])]
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
        assert qend>qstart
        if flag & 16:
            strand = '-'
        else:
            strand = '+'
        tname = line[2]
        tlen = lengths[tname]
        tstart = int(line[3])-1
        tend = tstart + 1

        fixed_cigar = list()
        fixed_cigar.append(cigar[0])
        for c,t in cigar[1:]:
            last_c,last_t = fixed_cigar[-1]
            if t==last_t:
                fixed_cigar[-1]=(last_c+c,t)
                continue
            if t in ['D','N'] and last_t in ['D','N']:
                fixed_cigar[-1]=(c+last_c,'N')
                continue
            fixed_cigar.append((c,t))
        cigar = fixed_cigar
        for c,t in cigar:
            if t in ['D','M','=','X']:
                tend+=c
        qstart_c = qstart
        qend_c   = qstart
        tstart_c = tstart
        tend_c   = tstart

        intervals = list()
        interval_cigar = list()
        for c,t in cigar:
            if t in ['D','I','X', '=', 'M']:
                interval_cigar.append('{}{}'.format(c,t))
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
                    ''.join(interval_cigar)
                ))
                interval_cigar = list()
                tend_c+=c
                tstart_c=tend_c
                qstart_c=qend_c
        if tstart_c < tend_c:
            intervals.append((
                qstart_c,
                qend_c,
                tstart_c,
                tend_c,
                ''.join(interval_cigar)
            ))
        print(qlen)
        for qstart_c,qend_c,tstart_c,tend_c,interval_cigar in intervals:
            print(qstart_c,qend_c,tstart_c,tend_c,)
            if strand == '-':
                qstart_c = qlen - qstart_c
                qend_c   = qlen - qend_c
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
                'oc:c:1',
                'cg:Z:{}'.format(interval_cigar),
            ]), file=outfile)
    outfile.close()

def main():
    args = parse_args()
    output_paf_from_sam(sam=args.sam, paf=args.paf)

if __name__ == "__main__":
    main()
