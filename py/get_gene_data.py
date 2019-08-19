#!/usr/bin/env python3
import argparse
import os
import sys
import pyfaidx
import pysam
import numpy as np
from Bio.Seq import Seq

POLY_A_LEN = 20
CIGAR = {
    0 : 'MATCH',
    1 : 'INS',
    2 : 'DEL',
    3 : 'REF_SKIP',
    4 : 'SOFT_CLIP',
    5 : 'HARD_CLIP',
    6 : 'PAD',
    7 : 'EQUAL',
    8 : 'DIFF',
    9 : 'BACK',
}
def parse_args():
    parser = argparse.ArgumentParser(
        description="Output the given gene, its reads, and transcripts. Any intronic region will be omitted")
    parser.add_argument("-t",
                        "--gtf",
                        type=str,
                        default='/groups/hachgrp/annotations/GTF/92/Homo_sapiens.GRCh38.92.gtf',
                        help="Path to GTF annotation")
    parser.add_argument("-d",
                        "--dna",
                        type=str,
                        default='/groups/hachgrp/annotations/DNA/92/Homo_sapiens.GRCh38.dna.fa',
                        help="Path to DNA genome reference")
    parser.add_argument("--minimap2",
                        type=str,
                        default='minimap2',
                        help="Path to minimap2")
    parser.add_argument("-r",
                        "--reads",
                        type=str,
                        required=True,
                        help="Path to reads in FASTQ format or in SAM alignment to genome reference")
    parser.add_argument("-g",
                        "--gene",
                        type=str,
                        required=True,
                        help="Gene name or ENSEMBL ID")
    parser.add_argument("-p",
                        "--padding",
                        type=int,
                        default=1000,
                        help="Pad size before and after the gene")
    parser.add_argument("-c",
                        "--threads",
                        type=int,
                        default=8,
                        help="Number of threads for minimap2. Default: 32")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        help="Output directory")
    parser.add_argument("-so",
                        "--sam-output",
                        type=str,
                        default="",
                        help="Output path for SAM file. Default is genome.sam under --output directory")
    args = parser.parse_args()
    return args

def query_yes_no(question, default=None):
    """Ask a yes/no question via input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".

    Copied from: https://stackoverflow.com/questions/3041986/apt-command-line-interface-like-yes-no-input
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

def get_coordinates_zero_based_end_exclusive(gene, gtf):
    if (gene[0:4]=='ENSG'):
        try:
            int(gene[4:])
            gene_gtf_id = 'gene_id "{}"'.format(gene)
        except:
            gene_gtf_id = 'gene_name "{}"'.format(gene)
    else:
        gene_gtf_id = 'gene_name "{}"'.format(gene)

    print('Getting {} coordinates from {}'.format(gene, gtf))
    gene_info = dict()
    found = False
    for line in open(gtf):
        if (line[0] == '#'):
            continue
        line = line.rstrip().split('\t')
        if (line[2] != 'gene'):
            continue
        if (gene_gtf_id in line[8]):
            if len(gene_info) == 0:
                info = {x.split()[0] : x.split()[1].strip('"') for x in line[8].strip('; ').split(';')}
                gene_info['gene_id']   = info['gene_id']
                gene_info['gene_name'] = info['gene_name']
                gene_info['chr']       = line[0]
                gene_info['strand']    = line[6]
                one_based_start_inc    = int(line[3])
                zero_based_start_inc = one_based_start_inc - 1
                gene_info['start'] = zero_based_start_inc
                one_based_end_inc      = int(line[4])
                zero_based_end_exl = one_based_end_inc
                gene_info['end'] = zero_based_end_exl
                found = True
                result = line
            else:
                print("There are multiple gene records in GTF with GTF gene identifier {}. We will pick only the first one. Offending line:".format(gene_gtf_id))
                print(line)
    if found:
        print(result)
        return gene_info
    else:
        print('Could not find GTF gene identifier {} in {} GTF file'.format(gene, gtf))
        exit(-1)

def align_reads(minimap2, threads, genome, reads, sam):
    if os.path.isfile(sam):
        if not query_yes_no(question='====\nSAM output file already exists "{}".\nDo you want to overwrite it?'.format(sam)):
            print('Existing..')
            exit()
    print('Aligning reads using command:')
    cmd = "{} -aY -x splice -t {} --secondary=no {} {} > {}".format(minimap2, threads, genome, reads, sam)
    print(cmd)
    os.system(cmd)

def get_exonic_positions(read, padding):
    i = read.reference_start
    for pos in range(i-padding,i):
        yield pos
    for (op,l) in read.cigartuples:
        op = CIGAR[op]
        if op in ['MATCH', 'DEL', 'EQUAL', 'DIFF']:
            for pos in range(i,i+l):
                yield pos
            i += l
        elif op in ['INS', 'SOFT_CLIP']:
            continue
        elif op in ['REF_SKIP']:
            i += l
            for pos in range(i-padding,i):
                yield pos
        else:
            raise Exception('Unrecognized cigar op: {}'.format(op))
    i = read.reference_end
    for pos in range(i,i+padding):
        yield pos

def get_gene_reads(sam, gene_info, filter_out_path):
    result = list()
    chr = gene_info['chr']
    strand = gene_info['strand']
    start = gene_info['padded_start']
    end   = gene_info['padded_end']
    print('Outputing {} reads of the gene at coordinates ({}) {}-{}'.format(sam, chr, start, end))

    sam = pysam.AlignmentFile(sam).fetch(contig=chr, start=start, stop=end)

    filter_out = open(filter_out_path, 'w+')
    for idx,read in enumerate(sam):
        to_filter_out = False
        if read.flag > 255:
            to_filter_out = True
        if read.mapping_quality < 60:
            to_filter_out = True
        if read.reference_name != chr:
            to_filter_out = True
        if read.reference_start < start:
            to_filter_out = True
        if read.reference_end > end:
            to_filter_out = True
        if to_filter_out:
            print('>{} {}:{}-{}'.format(read.qname, read.reference_name, read.reference_start, read.reference_end), file=filter_out)
            print(read.query, file=filter_out)
        else:
            result.append(read)
        if (idx+1)%100==0:
            print('Added {}/{} reads'.format(len(result),idx+1))
    filter_out.close()
    return result

def output_gene(genome, gene_info, out_path, polyA_len=POLY_A_LEN):
    gene_id              = gene_info['gene_id']
    gene_name            = gene_info['gene_name']
    chr                  = gene_info['chr']
    strand               = gene_info['strand']
    start                = gene_info['padded_start']
    end                  = gene_info['padded_end']
    is_exonic            = gene_info['is_exonic']
    genomic_to_genic_pos = gene_info['genomic_to_genic_pos']
    assert end-start == gene_info['padded_len']
    print('Outputting {} ({}) from reference {}'.format(gene_name, gene_id, genome))

    is_exonic_keys,is_exonic_vals = [list(l) for l in zip(*sorted(is_exonic.items()))]
    exon_idx_intervals = get_stretches_of_x(l=is_exonic_vals, x=True)
    print(exon_idx_intervals)

    genome = pyfaidx.Fasta(genome)
    out_file = open(out_path, 'w+')
    print('>{}.{}'.format(gene_id, gene_name), file=out_file)
    if (strand == '+'):
        for start,end in exon_idx_intervals:
            print((start, end), (is_exonic_keys[start], is_exonic_keys[end-1]))
            exon_seq = genome[chr][is_exonic_keys[start]:is_exonic_keys[end-1]]
            print(exon_seq, file=out_file, end='')
    elif (strand == '-'):
        for start,end in reversed(exon_idx_intervals):
            print((start, end), (is_exonic_keys[start], is_exonic_keys[end-1]))
            exon_seq = str(Seq(str(genome[chr][is_exonic_keys[start]:is_exonic_keys[end-1]])).reverse_complement())
            print(exon_seq, file=out_file, end='')
    else:
        raise Exception('strand_error:"{}"'.format(strand))
    print('A'*polyA_len, file=out_file)
    out_file.close()

def get_transcript_info(gtf, gene_info):
    gene_id = gene_info['gene_id']
    gene_name = gene_info['gene_name']
    chr = gene_info['chr']
    strand = gene_info['strand']
    print('Outputting {} transcripts from {}'.format(gene_name, gtf))

    transcript_infos = dict()
    for line in open(gtf):
        if (line[0] == '#'):
            continue
        line = line.rstrip().split('\t')
        if (line[2] != 'exon'):
            continue
        if (not gene_id in line[8]):
            continue
        info = {x.split()[0] : x.split()[1].strip('"') for x in line[8].strip('; ').split(';')}
        if (info['gene_biotype'] != 'protein_coding'):
            continue

        one_based_start_inc  = int(line[3])
        zero_based_start_inc = one_based_start_inc - 1
        exon_start = zero_based_start_inc

        one_based_end_inc = int(line[4])
        zero_based_end_exl = one_based_end_inc
        exon_end = zero_based_end_exl

        key = info['transcript_id']
        if key in transcript_infos:
            transcript_infos[key].append([exon_start, exon_end])
        else:
            transcript_infos[key] = [[exon_start, exon_end]]
    return transcript_infos

def output_transcripts(genome, transcript_infos, gene_info, out_tsv_path, out_seq_path):
    gene_id   = gene_info['gene_id']
    gene_name = gene_info['gene_name']
    chr       = gene_info['chr']
    strand    = gene_info['strand']
    genome    = pyfaidx.Fasta(genome)

    out_tsv = open(out_tsv_path, 'w+')
    out_seq = open(out_seq_path, 'w+')
    for tid, exons in transcript_infos.items():
        sorted(exons, reverse=strand=='-')
        print('{}\t{}\t{}\t'.format(tid, chr, strand), file=out_tsv, end='')
        print('>{}'.format(tid), file=out_seq)
        for start,end in exons:
            exon = str(genome[chr][start:end])
            assert(len(exon) == end-start)
            if (strand == '+'):
                interval = '{}-{}'.format(gene_info['genomic_to_genic_pos'][start], gene_info['genomic_to_genic_pos'][end])
            elif (strand == '-'):
                exon = str(Seq(exon).reverse_complement())
                interval = '{}-{}'.format(gene_info['genomic_to_genic_pos'][end], gene_info['genomic_to_genic_pos'][start])
            else:
                raise Exception('Unknown strand {}'.format(strand))
            print('{}'.format(exon), file=out_seq, end='')
            if [start,end]==exons[-1]:
                print('{}'.format(interval), file=out_tsv, end='')
            else:
                print('{}'.format(interval), file=out_tsv, end=',')
        print('', file=out_tsv)
        print('', file=out_seq)
    out_tsv.close()
    out_seq.close()

def get_chr_lens(fai):
    result = dict()
    for line in open(fai):
        line = line.rstrip().split('\t')
        assert len(line) == 5, 'FAI records each has 5 fields'
        result[line[0]] = int(line[1])
    return result

def get_stretches_of_x(l, x):
    result = list()
    idx = 0
    while idx < len(l):
        s = idx
        while s < len(l) and l[s] != x:
            s+=1
        e = s + 1
        while e < len(l) and l[e] == x:
            e+=1
        if e < len(l) and l[e] == x:
            e+=1
        result.append((s,e))
        idx = e + 1
    return result

def output_reads(reads, gene_info, out_path):
    out_file = open(out_path, 'w+')
    for read in reads:
        print('>{} {}:{}-{}'.format(read.qname, read.reference_name, read.reference_start, read.reference_end), file=out_file)
        if (gene_info['strand'] == '+'):
            print(read.query_sequence, file=out_file)
        elif (gene_info['strand'] == '-'):
            print(Seq(read.query_sequence).reverse_complement(), file=out_file)
        else:
            raise Exception('Unknown strand {}'.format(gene_info['strand']))

def main():
    args = parse_args()
    args.output = args.output.rstrip('/')
    args.output += '/'
    print('Making output directory: {}'.format(args.output))
    os.makedirs(args.output, exist_ok=True)

    gene_info = get_coordinates_zero_based_end_exclusive(gene=args.gene, gtf=args.gtf)
    chr_lens = get_chr_lens(fai='{}.fai'.format(args.dna))
    gene_info['padded_start'] = max(gene_info['start'] - args.padding, 0)
    gene_info['padded_end']   = min(gene_info['end']   + args.padding, chr_lens[gene_info['chr']])
    gene_info['padded_len'] = gene_info['padded_end'] - gene_info['padded_start']

    is_exonic = {pos:False for pos in range(gene_info['padded_start'],gene_info['padded_end'])}
    if (args.reads):
        try:
            sam = pysam.AlignmentFile(args.reads)
            sam.close()
            print('Reads {} file is in SAM/BAM format'.format(args.reads))
            sam = args.reads
        except ValueError:
            print('Reads {} file is NOT in SAM format. Assuming FASTA/Q format'.format(args.reads))
            if (args.sam_output == ""):
                sam = '{}genome.sam'.format(args.output)
            else:
                sam = args.sam_output
            align_reads(minimap2=args.minimap2, threads=args.threads, genome=args.dna, reads=args.reads, sam=sam)
        reads = get_gene_reads(sam=sam, gene_info=gene_info, filter_out_path='{}reads.filtered_out.fasta'.format(args.output))
        print('There are {} reads that belong to the gene'.format(len(reads)))
        output_reads(reads=reads, gene_info=gene_info, out_path='{}reads.fasta'.format(args.output))
        for read in reads:
            for pos in get_exonic_positions(read=read, padding=args.padding):
                is_exonic[pos] = True

    transcript_infos = get_transcript_info(gtf=args.gtf, gene_info=gene_info)
    for tid,exons in transcript_infos.items():
        for (start,end) in exons:
            for pos in range(start-args.padding,end+args.padding):
                is_exonic[pos] = True

    genomic_positions = sorted(is_exonic.keys(), reverse=gene_info['strand']=='-')
    genomic_to_genic_pos = dict()
    genic_pos = 0
    for pos in genomic_positions:
        genic_pos+=is_exonic[pos]
        genomic_to_genic_pos[pos] = genic_pos
    gene_info['genomic_to_genic_pos'] = genomic_to_genic_pos
    gene_info['is_exonic'] = is_exonic

    output_transcripts(
        genome           = args.dna,
        transcript_infos = transcript_infos,
        gene_info        = gene_info,
        out_tsv_path     = '{}transcripts.tsv'.format(args.output),
        out_seq_path     = '{}transcripts.fasta'.format(args.output)
    )
    output_gene(
        genome    = args.dna,
        gene_info = gene_info,
        out_path  = '{}gene.fasta'.format(args.output)
    )

if __name__ == "__main__":
    main()
