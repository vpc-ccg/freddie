#!/usr/bin/env python3
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster all barcodes with hamming distance thershold")
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
                        "--pad-size",
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

import sys

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
    print('Getting {} coordinates from {}'.format(gene, gtf))
    gene_info = dict()
    found = False
    for line in open(gtf):
        if (line[0] == '#'):
            continue
        line = line.rstrip().split('\t')
        if (line[2] != 'gene'):
            continue
        if (gene in line[8]):
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
                print("There are multiple gene records in GTF with GTF gene identifier {}. We will pick only the first one. Offending line:".format(gene))
                print(line)
    if found:
        print(result)
        return gene_info
    else:
        print('Could not find GTF gene identifier {} in {} GTF file'.format(gene, gtf))
        exit(-1)

def align_reads(minimap2, threads, genome, reads, sam):
    import os
    if os.path.isfile(sam):
        if not query_yes_no(question='====\nSAM output file already exists "{}".\nDo you want to overwrite it?'.format(sam)):
            print('Existing..')
            exit()
    print('Aligning reads using command:')
    cmd = "{} -aY -x splice -t {} --secondary=no {} {} > {}".format(minimap2, threads, genome, reads, sam)
    print(cmd)
    import os
    os.system(cmd)

def output_gene_reads(sam, gene_info, padding, out):
    chr = gene_info['chr']
    strand = gene_info['strand']
    start = gene_info['start']
    end = gene_info['end']
    print('Outputing {} reads of the gene at coordinates ({}) {}-{} and left padding of {}'.format(sam, chr, start, end, padding))
    import pysam
    from Bio.Seq import Seq
    sam = pysam.AlignmentFile(sam)
    if (sam.has_index()):
        print('{} has an index. Using that'.format(sam))
        sam = sam.fetch(contig=chr, start=max(start-padding,0), stop=end)
    else:
        print('{} has no index. Scanning the whole thing'.format(sam))
        sam = sam.fetch(until_eof=True)

    print('>{} {}:{}-{}'.format('UNALINED_READ', 'NO_CHROM', '-1', '-1'), file=out)
    print('NNNNNNNNNNNNNNNNNNNNNNN', file=out)
    for idx,read in enumerate(sam):
        if idx%500==0:
            print('Read the {}-th read'.format(idx))
        if (read.flag > 255 or read.mapping_quality < 60):
            continue
        if (read.reference_name == chr and read.reference_start > start-padding and read.reference_start < end):
            print('>{} {}:{}-{}'.format(read.qname, read.reference_name, read.reference_start, read.reference_end), file=out)
            if (strand == '+'):
                print(read.query, file=out)
            elif (strand == '-'):
                print(Seq(read.query).reverse_complement(), file=out)
            else:
                print('strand_error:"{}"'.format(strand), file=out)

def output_gene(genome, gene_info, out):
    gene_id = gene_info['gene_id']
    gene_name = gene_info['gene_name']
    chr = gene_info['chr']
    strand = gene_info['strand']
    start = gene_info['start']
    end = gene_info['end']
    print('Outputting {} ({}) from reference {}'.format(gene_name, gene_id, genome))
    import pyfaidx
    from Bio.Seq import Seq
    genome = pyfaidx.Fasta(genome)
    print('>{}.{}'.format(gene_id, gene_name), file=out)
    if (strand == '+'):
        print(genome[chr][start:end], file=out)
    elif (strand == '-'):
        print(str(Seq(str(genome[chr][start:end])).reverse_complement()), file=out)
    else:
        print('strand_error:"{}"'.format(strand), file=out)

def output_transcripts(genome, gtf, gene_info, out_tsv, out_seq):
    gene_id = gene_info['gene_id']
    gene_name = gene_info['gene_name']
    chr = gene_info['chr']
    strand = gene_info['strand']
    gene_start = gene_info['start']
    gene_end = gene_info['end']
    print('Outputting {} transcripts from {}'.format(gene_name, gtf))
    import pyfaidx
    from Bio.Seq import Seq

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

    genome = pyfaidx.Fasta(genome)
    for tid, exons in transcript_infos.items():
        if (strand == '+'):
            sorted(exons)
        else:
            sorted(exons, reverse=True)
        print('{}\t{}\t{}\t'.format(tid, chr, strand), file=out_tsv, end='')
        print('>{}'.format(tid), file=out_seq)
        for start,end in exons:
            exon = str(genome[chr][start:end])
            assert(len(exon) == end-start)
            if (strand == '+'):
                interval = '{}-{}'.format(start-gene_start, end-gene_start)
            elif (strand == '-'):
                exon = str(Seq(exon).reverse_complement())
                interval = '{}-{}'.format(gene_end-end, gene_end-start)
            else:
                exon = 'err_strand'
                interval = strand
            print('{}'.format(exon), file=out_seq, end='')
            if [start,end]==exons[-1]:
                print('{}'.format(interval), file=out_tsv, end='')
            else:
                print('{}'.format(interval), file=out_tsv, end=',')
        print('', file=out_tsv)
        print('', file=out_seq)

def main():
    args = parse_args()
    import os
    args.output = args.output.rstrip('/')
    args.output += '/'
    print('Making output directory: {}'.format(args.output))
    os.makedirs(args.output, exist_ok=True)
    if (args.gene[0:4]=='ENSG'):
        try:
            int(args.gene[4:])
            gene_gtf_id = 'gene_id "{}"'.format(args.gene)
        except:
            gene_gtf_id = 'gene_name "{}"'.format(args.gene)
    else:
        gene_gtf_id = 'gene_name "{}"'.format(args.gene)

    gene_info = get_coordinates_zero_based_end_exclusive(gene=gene_gtf_id, gtf=args.gtf)

    if (args.reads):
        import pysam
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
        out = open('{}reads.fasta'.format(args.output), 'w+')
        output_gene_reads(sam=sam, gene_info=gene_info, padding=args.pad_size, out=out)
        out.close()

    out = open('{}gene.fasta'.format(args.output), 'w+')
    output_gene(genome=args.dna, gene_info=gene_info, out=out)
    out.close()

    out_tsv = open('{}transcripts.tsv'.format(args.output), 'w+')
    out_seq = open('{}transcripts.fasta'.format(args.output), 'w+')
    output_transcripts(genome=args.dna, gtf=args.gtf, gene_info=gene_info, out_tsv=out_tsv, out_seq=out_seq)
    out_tsv.close()
    out_seq.close()

if __name__ == "__main__":
    main()
