#!/usr/bin/env python3
import sys

paf = sys.argv[1]
target = sys.argv[2]
output = sys.argv[3]
contig = sys.argv[4]
tid = sys.argv[5]
 
iso_seq = open(target).readlines()[1].rstrip()
locs = [0 for _ in range(len(iso_seq))]
rids = set()
for l in open(paf):
    l = l.rstrip().split('\t')
    rids.add(l[0])
    assert len(locs) == int(l[6])
    for i in range(int(l[7]),int(l[8])):
        locs[i]+=1

out_file = open(output, 'w+')
out_file.write('>isoform_{}_{}\n'.format(contig, tid))
for i,c in zip(locs,iso_seq):
    if len(rids) == 0:
        break
    if i/len(rids) > 0.3:
        out_file.write(c)
out_file.write('\n')
out_file.close()
