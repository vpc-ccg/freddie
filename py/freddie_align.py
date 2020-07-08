#!/usr/bin/env python3
import argparse
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(
        description="Outputs SAM alignments using deSALT splice aligner")
    parser.add_argument("-r",
                        "--reads",
                        nargs="+",
                        type=str,
                        required=True,
                        help="Space separated paths to reads in FASTQ or FASTA format")
    parser.add_argument("-g",
                        "--genome",
                        type=str,
                        default=None,
                        help="Path to FASTA genome reference file. Needed if --desalt-index is not provided.")
    parser.add_argument("-od",
                        "--out-desalt-index",
                        type=str,
                        default=None,
                        help="Path to output deSALT index directory. Needed if --desalt-index is not provided.")
    parser.add_argument("-d",
                        "--desalt",
                        type=str,
                        default='deSALT',
                        help="Path to deSALT executable. Default: deSALT")
    parser.add_argument("-i",
                        "--desalt-index",
                        type=str,
                        default=None,
                        help="Path to deSALT index directory. If provided, --genome and --out-desalt-index will not be used")
    parser.add_argument("-t",
                        "--threads",
                        type=int,
                        default=1,
                        help="Number of threads for deSALT alignment. Default: 1")
    parser.add_argument("-m",
                        "--temporary-prefix",
                        type=str,
                        default=None,
                        help="Prefix for deSALT temporary alignment files. Default: <--output>.temp")
    parser.add_argument("-s",
                        "--sequencer",
                        type=str,
                        default='null',
                        help="Type of sequencing machine to provide deSALT: null, ccs, clr, ont1d, or ont2d. Default: null")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        help="Output file for deSALT SAM alignments")
    args = parser.parse_args()
    assert args.sequencer in ['null', 'ccs', 'clr', 'ont1d', 'ont2d']
    if args.desalt_index == None:
        assert args.genome != None, '--genome must be provided since --desalt-index is missing!'
        assert args.out_desalt_index != None, '--out-desalt-index must be provided since --desalt-index is missing!'
    if args.temporary_prefix == None:
        args.temporary_prefix = '{}.temp'.format(args.output)

    return args

def main():
    args = parse_args()
    subprocess.run(['pwd'])

    if args.desalt_index == None:
        args.desalt_index = args.out_desalt_index
        run_args = list()
        run_args.append(args.desalt)
        run_args.append('index')
        run_args.append(args.genome)
        run_args.append(args.desalt_index)
        print('Running: {}'.format(' '.join(run_args)))
        completed = subprocess.run(run_args)
        assert completed.returncode == 0

    run_args = list()
    run_args.append(args.desalt)
    run_args.append('aln')
    run_args.append('-f')
    run_args.append(args.temporary_prefix)
    if args.sequencer != 'null':
        run_args.append('-x')
        run_args.append(args.sequencer)
    run_args.append('-t')
    run_args.append(str(args.threads))
    run_args.append('-o')
    run_args.append(args.output)
    run_args.append(args.desalt_index)
    run_args.extend(args.reads)
    print('Running: {}'.format(' '.join(run_args)))
    completed = subprocess.run(run_args)
    assert completed.returncode == 0

if __name__ == "__main__":
    main()
