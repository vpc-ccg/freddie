#!/usr/bin/env python
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
                        default=None,
                        help="Path to reads in FASTQ format or in SAM alignment to genome reference")
    parser.add_argument("-g",
                        "--gene",
                        type=str,
                        required=True,
                        help="Gene name or ENSEMBL ID")
    parser.add_argument("-p",
                        "--pad-size",
                        type=int,
                        default=100,
                        help="Pad size before and after the gene")
    parser.add_argument("-temp",
                        "--threads",
                        type=int,
                        default=32,
                        help="Number of threads for minimap2. Default: 32")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        help="Output directory")
    args = parser.parse_args()
    return args

def get_coordinates(gene, gtf):
    print('Getting {} coordinates from {}'.format(gene, gtf))
    for line in open(gtf):
        if (line[0] == '#'):
            continue
        line = line.rstrip().split('\t')
        if (line[2] != 'gene'):
            continue
        if (gene in line[8]):
            print(line)
            return line[0], int(line[3]), int(line[4])

def align_reads(minimap2, threads, genome, reads, sam):
    import os
    os.system("{} -aY -x splice -k 14 -t {} --secondary=no {} {} > {}".format(minimap2, threads, genome, reads, sam))

def output_gene_reads(sam, chr, start, end, padding, out):
    import pysam
    sam = pysam.AlignmentFile(sam)
    if (sam.has_index()):
        sam = sam.fetch(contig=chr, start=max(start-padding,0), stop=end)
    else:
        sam = sam.fetch(until_eof=True)
    for read in sam:
        print('>{}\n{}'.format(read.qname, read.query), file=out)

def output_gene(genome, gene, chr, start, end, out):
    print('Outputting {} from reference {}'.format(gene, genome))
    import pyfaidx
    from Bio.Seq import Seq
    genome = pyfaidx.Fasta(genome)
    print('>{}\n{}'.format(gene, genome[chr][start:end]), file=out)

def output_transcripts(genome, gtf, gene, gene_start, out_tsv, out_seq):
    print('Outputting {} transcripts from {}'.format(gene, gtf))
    import pyfaidx
    from Bio.Seq import Seq

    transcript_infos = dict()
    for line in open(gtf):
        if (line[0] == '#'):
            continue
        line = line.rstrip().split('\t')
        if (line[2] != 'exon'):
            continue
        if (not gene in line[8]):
            continue
        info = {x.split()[0] : x.split()[1].strip('"') for x in line[8].strip('; ').split(';')}
        if (info['gene_biotype'] != 'protein_coding'):
            continue

        exon_start, exon_end = int(line[3]), int(line[4])
        chr, strand = line[0], line[6]
        key = (info['transcript_id'], chr, strand)
        if key in transcript_infos:
            transcript_infos[key].append([exon_start, exon_end])
        else:
            transcript_infos[key] = [[exon_start, exon_end]]

    genome = pyfaidx.Fasta(genome)
    for (tid, chr, strand), exons in transcript_infos.items():
        if (strand == '+'):
            sorted(exons)
        else:
            sorted(exons, reverse=True)
        print('{}\t{}\t{}\t'.format(tid, chr, strand), file=out_tsv, end='')
        print('>{}'.format(tid), file=out_seq)
        for start,end in exons:
            if (strand == '+'):
                exon = str(Seq(str(genome[chr][start:end])))
                interval = '{}-{}'.format(start-gene_start, end-gene_start)
            else:
                exon = str(Seq(str(genome[chr][start:end])).reverse_complement())
                interval = '{}-{}'.format(end-gene_start, start-gene_start)
            print('{}'.format(exon), file=out_seq, end='')
            print('{}'.format(interval), file=out_tsv, end=',')
        print('', file=out_tsv)
        print('', file=out_seq)

def main():
    args = parse_args()
    if (args.gene[0:4]=='ENSG'):
        try:
            int(args.gene[4:])
            args.gene = 'gene_id "{}"'.format(args.gene)
        except:
            args.gene = 'gene_name "{}"'.format(args.gene)
    else:
        args.gene = 'gene_name "{}"'.format(args.gene)

    chr, start, end = get_coordinates(args.gene, args.gtf)

    if (args.reads):
        import pysam
        try:
            sam = pysam.AlignmentFile(args.reads)
            sam.close()
            sam = args.alignments
        except ValueError:
            sam = '{}genome.sam'.format(args.output)
            align_reads(args.minimap2, args.dna, args.threads, args.reads, sam)
        out = open('{}reads.fasta'.format(args.output), 'w+')
        output_gene_reads(sam, chr, start, end, args.padding, out)
        out.close()

    out = open('{}gene.fasta'.format(args.output), 'w+')
    output_gene(args.dna, args.gene, chr, start, end, out)
    out.close()

    out_tsv = open('{}trasncripts.tsv'.format(args.output), 'w+')
    out_seq = open('{}trasncripts.fasta'.format(args.output), 'w+')
    output_transcripts(args.dna, args.gtf, args.gene, start, out_tsv, out_seq)
    out_tsv.close()
    out_seq.close()

if __name__ == "__main__":
    main()
