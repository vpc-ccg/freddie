#!/usr/bin/env python3
import argparse
import pyfaidx
import os

def parse_args():
    parser = argparse.ArgumentParser(
        description="Output the given gene, its reads, and transcripts. Any intronic region will be omitted")
    parser.add_argument("-a",
                        "--genome",
                        type=str,
                        default='/groups/hachgrp/annotations/DNA/97/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
                        help="Path to DNA genome reference")
    parser.add_argument("-t",
                        "--gtf",
                        type=str,
                        default='/groups/hachgrp/annotations/GTF/97/Homo_sapiens.GRCh38.97.gtf',
                        help="Path to GTF annotation")
    parser.add_argument("-q",
                        "--fastqs",
                        nargs="+",
                        type=str,
                        required=True,
                        help="Space separated paths to reads in FASTQ format")
    parser.add_argument("-g",
                        "--gene",
                        type=str,
                        required=True,
                        help="Gene name or ENSEMBL ID")
    parser.add_argument("-p",
                        "--paf",
                        type=str,
                        required=True,
                        help="Path to PAF file of aligned reads")
    parser.add_argument("-d",
                        "--padding",
                        type=int,
                        default=1000,
                        help="Pad size before and after the gene")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        help="Output directory")
    args = parser.parse_args()
    return args

def get_transcript_info(gtf, gene_id):
    print('Outputting {} transcripts from {}'.format(gene_id, gtf))
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
        # if (info['gene_biotype'] != 'protein_coding'):
        #     continue
        one_based_start_inc  = int(line[3])
        zero_based_start_inc = one_based_start_inc - 1
        exon_start = zero_based_start_inc

        one_based_end_inc = int(line[4])
        zero_based_end_exl = one_based_end_inc
        exon_end = zero_based_end_exl

        key = info['transcript_id']
        if not key in transcript_infos:
            transcript_infos[key] = list()
        transcript_infos[key].append((exon_start, exon_end))
    return transcript_infos

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
    found = 0
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

def output_gene(genome_fasta, gene_info, out_path):
    genome = pyfaidx.Fasta(genome_fasta)
    out_file = open('{}gene.fasta'.format(out_path, ), 'w+')
    print('>{gene_id}.{gene_name} chr={chr} final_start={final_start} final_end={final_end} gene_start={gene_start} gene_end={gene_end}'.format(
        gene_id=gene_info['gene_id'],
        gene_name=gene_info['gene_name'],
        chr=gene_info['chr'],
        final_start=gene_info['final_start'],
        final_end=gene_info['final_end'],
        gene_start=gene_info['start']-gene_info['final_start'],
        gene_end=gene_info['end']-gene_info['final_start'],
    ), file=out_file)
    print(genome[gene_info['chr']][gene_info['final_start']:gene_info['final_end']], file=out_file)
    out_file.close()

def get_gene_reads(gene_info, paf, fastqs, out_path):
    rids = set()
    start = gene_info['padded_start']
    end   = gene_info['padded_end']
    for line in open(paf):
        fields = line.rstrip().split('\t')
        tstart = int(fields[7])
        tend   = int(fields[8])
        if fields[5]==gene_info['chr'] and (tstart <= gene_info['padded_start'] <= tend or gene_info['padded_start'] <= tstart <= gene_info['padded_end']):
            rids.add(fields[0])
    for line in open(paf):
        fields = line.rstrip().split('\t')
        if not fields[0] in rids:
            continue
        start = min(start, int(fields[7]))
        end   = max(end, int(fields[8]))
    out_paf = open('{}reads.paf'.format(out_path), 'w+')
    for line in open(paf):
        fields = line.rstrip().split('\t')
        if not fields[0] in rids:
            continue
        fields[5] = str(gene_info['gene_name'])
        fields[6] = str(end - start)
        fields[7] = str(int(fields[7]) - start)
        fields[8] = str(int(fields[8]) - start)
        print('\t'.join(fields), file=out_paf)
    out_paf.close()
    out_fastq = open('{}reads.fastq'.format(out_path),'w+')
    for fq in fastqs:
        for line in open(fq):
            if line[0]=='@':
                if line.rstrip().split()[0][1:] in rids:
                    print_flag = True
                else:
                    print_flag = False
            if print_flag:
                out_fastq.write(line)
    out_fastq.close()
    gene_info['final_start'] = start
    gene_info['final_end']   = end

def output_transcripts(gene_info, transcript_infos, genome_fasta, out_path):
    out_file = open('{}transcripts.tsv'.format(out_path),'w+')
    start = gene_info['final_start']
    for tid,genome_intervals in transcript_infos.items():
        print('\t'.join([
            tid,gene_info['chr'],
            gene_info['strand'],
            ','.join(['{}-{}'.format(gstart-start,gend-start) for gstart,gend in genome_intervals]),
            ','.join(['{}-{}'.format(gstart,gend) for gstart,gend in genome_intervals]),
        ]), file=out_file)
    out_file.close()

    out_file = open('{}transcripts.fasta'.format(out_path),'w+')
    genome = pyfaidx.Fasta(genome_fasta)
    for tid,genome_intervals in transcript_infos.items():
        print('>{} genome_intervals={}'.format(
            tid,
            ','.join(['{}-{}'.format(gstart,gend) for gstart,gend in genome_intervals]),
        ),file=out_file)
        print(''.join([str(genome[gene_info['chr']][gstart:gend]) for gstart,gend in genome_intervals]), file=out_file)
    out_file.close()

def main():
    args = parse_args()
    args.output = args.output.rstrip('/')
    args.output += '/'
    print('Making output directory: {}'.format(args.output))
    os.makedirs(args.output, exist_ok=True)

    gene_info = get_coordinates_zero_based_end_exclusive(gene=args.gene, gtf=args.gtf)
    gene_info['padded_start'] = max(gene_info['start'] - args.padding,0)
    gene_info['padded_end']   = gene_info['end'] + args.padding
    gene_info['padded_len'] = gene_info['padded_end'] - gene_info['padded_start']

    get_gene_reads(gene_info=gene_info, paf=args.paf, fastqs=args.fastqs, out_path=args.output)
    output_gene(genome_fasta=args.genome, gene_info=gene_info, out_path=args.output)
    transcript_infos = get_transcript_info(gtf=args.gtf, gene_id=gene_info['gene_id'])
    output_transcripts(gene_info=gene_info, transcript_infos=transcript_infos, genome_fasta=args.genome, out_path=args.output)

if __name__ == "__main__":
    main()
