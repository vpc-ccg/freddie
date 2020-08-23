#!/usr/bin/env python3
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(
        description="Converts GTF to BED files")
    parser.add_argument("-g",
                        "--gtf",
                        type=str,
                        required=True,
                        help="Path to GTF file of the transcripts")
    parser.add_argument("-c",
                        "--coverage-tsv",
                        type=str,
                        required=False,
                        help="Path to Freddie transcript coverage TSV file. If not provided, all transcripts will be output.")
    parser.add_argument("-od",
                        "--output-dir",
                        type=str,
                        default='freddie_beds/',
                        help="Path to output directory. Default: freddie_beds/")
    args = parser.parse_args()
    args.output_dir = args.output_dir.rstrip('/')
    return args

def main():
    args = parse_args()

    if args.coverage_tsv != None:
        tids = {line.split('\t')[0] for line in open(args.coverage_tsv)}

    os.makedirs(args.output_dir, exist_ok=True)
    for line in open(args.gtf):
        if line[0]=='#':
            continue
        line = line.rstrip().split('\t')
        if not line[2] in ['exon', 'transcript']:
            continue
        tid = [x.strip(' "').split('"')[1] for x in line[8].strip(';').split(';') if x.strip(' ').startswith('transcript_id')][0]
        if args.coverage_tsv != None and not tid in tids:
            continue
        if line[2] == 'transcript':
            if 'out_file' in locals():
                out_file.close()
            out_file = open('{}/{}.bed'.format(args.output_dir, tid), 'w+')
        if line[2] == 'exon':
            print(line[0], int(line[3])-1, line[4], tid, sep='\t', file=out_file)
    out_file.close()
if __name__ == "__main__":
    main()
