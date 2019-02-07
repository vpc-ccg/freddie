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
    parser.add_argument("-c",
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
            return line[0], line[6], int(line[3]), int(line[4])

def align_reads(minimap2, threads, genome, reads, sam):
    print('Aligning reads using command:')
    cmd = "{} -aY -x splice -t {} --secondary=no {} {} > {}".format(minimap2, threads, genome, reads, sam)
    print(cmd)
    import os
    os.system(cmd)

def output_gene_reads(sam, chr, strand, start, end, padding, out):
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
    for read in sam:
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

def output_gene(genome, gene, chr, strand, start, end, out):
    print('Outputting {} from reference {}'.format(gene, genome))
    import pyfaidx
    from Bio.Seq import Seq
    genome = pyfaidx.Fasta(genome)
    print('>{}'.format(gene), file=out)
    if (strand == '+'):
        print(genome[chr][start:end], file=out)
    elif (strand == '-'):
        print(str(Seq(str(genome[chr][start:end])).reverse_complement()), file=out)
    else:
        print('strand_error:"{}"'.format(strand), file=out)

def output_transcripts(genome, gtf, gene, chr, strand, gene_start, gene_end, out_tsv, out_seq):
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
            if (strand == '+'):
                exon = str(genome[chr][start:end])
                interval = '{}-{}'.format(start-gene_start, end-gene_start)
            elif (strand == '-'):
                exon = str(Seq(str(genome[chr][start:end])).reverse_complement())
                interval = '{}-{}'.format(gene_end-end, gene_end-start)
            else:
                exon = 'err_strand'
                interval = strand
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

    chr, strand, start, end = get_coordinates(gene=args.gene, gtf=args.gtf)

    if (args.reads):
        import pysam
        try:
            sam = pysam.AlignmentFile(args.reads)
            sam.close()
            print('Reads {} file is in SAM format'.format(args.reads))
            sam = args.reads
        except ValueError:
            print('Reads {} file is NOT in SAM format. Assuming FASTA/Q format'.format(args.reads))
            sam = '{}genome.sam'.format(args.output)
            align_reads(minimap2=args.minimap2, threads=args.threads, genome=args.dna, reads=args.reads, sam=sam)
        out = open('{}reads.fasta'.format(args.output), 'w+')
        output_gene_reads(sam=sam, chr=chr, strand=strand, start=start, end=end, padding=args.pad_size, out=out)
        out.close()

    out = open('{}gene.fasta'.format(args.output), 'w+')
    output_gene(genome=args.dna, gene=args.gene, chr=chr, strand=strand, start=start, end=end, out=out)
    out.close()

    out_tsv = open('{}trasncripts.tsv'.format(args.output), 'w+')
    out_seq = open('{}trasncripts.fasta'.format(args.output), 'w+')
    output_transcripts(genome=args.dna, gtf=args.gtf, gene=args.gene, chr=chr, strand=strand,  gene_start=start, gene_end=end, out_tsv=out_tsv, out_seq=out_seq)
    out_tsv.close()
    out_seq.close()

if __name__ == "__main__":
    main()
