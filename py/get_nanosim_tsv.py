#!/usr/bin/env python
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster all barcodes with hamming distance thershold")
    parser.add_argument("-nsr",
                        "--nanosim-reads",
                        type=str,
                        required=True,
                        help="Path to NanoSim reads in FASTQ format")
    parser.add_argument("-t",
                        "--transcript-tsv",
                        type=str,
                        required=True,
                        help="Transcripts TSV file")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        default=None,
                        help="Output file. Default stdout")
    args = parser.parse_args()
    return args

def get_transcript_infos(transcript_tsv):
    transcript_infos = dict()
    for line in open(transcript_tsv):
        line = line.rstrip().split('\t')
        tid = line[0]
        transcript_infos[tid] = dict()
        transcript_infos[tid]['chr']    = line[1]
        transcript_infos[tid]['strand'] = line[2]
        transcript_infos[tid]['exons'] = [(int(x.split('-')[0]),int(x.split('-')[1])) for x in line[3].rstrip(',').split(',')]
        transcript_infos[tid]['length'] = sum([x[1]-x[0] for x in transcript_infos[tid]['exons']])
    return transcript_infos

def output_nanosim_reads_tsv(transcript_infos, nanosim_reads_fasta, out_path):
    if out_path == None:
        import sys
        out_file = sys.stdout
    else:
        out_file = open(out_path, 'w+')
    for line in open(nanosim_reads_fasta):
        if line[0] != '>':
            continue
        rname = line[1:].rstrip().split()[0]
        comment = line[1:].rstrip().split()[1]
        line = rname.split('_')
        tid = line[0]
        start_on_transcript = int(line[1])
        aligned = line[2] == 'aligned'
        rid = int(line[3])
        strand = line[4]
        ssc = int(line[5])
        aln = int(line[6])
        esc = int(line[7])
        gene_starts = list()
        gene_ends = list()
        r_measured = 0
        t_measured = 0
        eidx = 0
        # Finding the first exon start
        for exon in transcript_infos[tid]['exons']:
            eidx += 1
            exon_len = exon[1]-exon[0]
            t_measured += exon_len
            if start_on_transcript < t_measured:
                offset = start_on_transcript - (t_measured - exon_len)
                gene_starts.append(exon[0] + offset)
                if aln - r_measured >= t_measured - start_on_transcript:
                    gene_ends.append(exon[1])
                    r_measured += t_measured - start_on_transcript
                else:
                    gene_ends.append(gene_starts[-1]+ (aln-r_measured))
                    r_measured += gene_ends[-1]-gene_starts[-1]
                break
        # From first exon start till last exon end
        for exon in transcript_infos[tid]['exons'][eidx:]:
            eidx+=1
            exon_len = exon[1]-exon[0]
            t_measured += exon_len
            if aln-r_measured <= 0:
                break
            gene_starts.append(exon[0])
            if aln-r_measured > exon_len:
                gene_ends.append(exon[1])
                r_measured+=exon_len
            else:
                gene_ends.append(gene_starts[-1] + (aln-r_measured))
                r_measured += gene_ends[-1]-gene_starts[-1]
        assert(r_measured == aln)
        for exon in transcript_infos[tid]['exons'][eidx:]:
            exon_len = exon[1]-exon[0]
            t_measured += exon_len
        assert(t_measured == transcript_infos[tid]['length'])
        exons = ','.join(('{}-{}'.format(x[0],x[1]) for x in zip(gene_starts, gene_ends)))
        print('\t'.join((rname, tid, strand,exons)), file=out_file)

def main():
    args = parse_args()
    transcript_infos = get_transcript_infos(args.transcript_tsv)
    output_nanosim_reads_tsv(transcript_infos, args.nanosim_reads, args.output)
if __name__ == "__main__":
    main()
