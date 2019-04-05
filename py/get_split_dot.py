#!/usr/bin/env python3
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Split each read annot in DOT file into a separate DOT file")
    parser.add_argument("-d",
                        "--dot",
                        type=str,
                        required=True,
                        help="Input DOT file")
    parser.add_argument("-t",
                        "--reads-tsv",
                        type=str,
                        required=True,
                        help="Transcripts TSV file")
    parser.add_argument("-o",
                        "--output-prefix",
                        type=str,
                        required=True,
                        help="Output prefix. Out files will be <prefix><read name>.dot")
    args = parser.parse_args()
    return args


def get_rname_to_tname(tsv_path):
    rname_to_tname = dict()
    for line in open(tsv_path):
        line = line.rstrip().split('\t')
        rname_to_tname[line[0]] = line[1]
    return rname_to_tname

def print_dot_per_read(rname_to_tname, dot_path, output_prefix):
    for rname in rname_to_tname:
        out_file = open('{}{}.dot'.format(output_prefix, rname), 'w+')
        for line in open(dot_path):
            line = line.rstrip();
            if not 'ANNOTATION' in line:
                print(line, file=out_file)
                continue
            seq_name = line.split()[-1]
            if rname == seq_name or rname_to_tname[rname] == seq_name:
                print(line, file=out_file)
                continue
        out_file.close()


def main():
    args = parse_args()
    rname_to_tname = get_rname_to_tname(args.reads_tsv)
    print_dot_per_read(rname_to_tname, args.dot, args.output_prefix)




if __name__ == "__main__":
    main()
