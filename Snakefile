configfile: 'config.yaml'

def get_abs_path(path):
    import os
    abs_path = os.popen('readlink -f {}'.format(path)).read()
    return abs_path.rstrip("\r\n")
def make_slurm():
    import os
    os.makedirs('slurm'.format(outpath), mode=0o777, exist_ok=True)

outpath = get_abs_path(config['outpath'])
make_slurm()
config['exec']['freddie'] = get_abs_path(config['exec']['freddie'])

genes_d    = '{}/genes'.format(outpath)
mapped_d   = '{}/mapped'.format(outpath)
ns_train_d = '{}/ns_train'.format(outpath)

gene_data=[
    'transcripts.tsv',
    'transcripts.fasta',
    'gene.fasta',
    'reads_raw.fasta',
    'reads_raw.filtered_out.fasta',
]
# ns_train_files=[
#     # '_unaligned_length.pkl',
#     # '_aligned_region.pkl',
#     # '_aligned_reads.pkl',
#     # '_ht_length.pkl',
#     # '_match.hist',
#     # '_mis.hist',
#     # '_del.hist',
#     # '_ins.hist',
#     # '_first_match.hist',
#     # '_error_markov_model',
#     # '_ht_ratio.pkl',
#     # '.sam',
#     # '_match_markov_model',
#     # '_model_profile',
#     # '_error_rate.tsv',
#     '_done'
# ]
#
rule all:
    input:
        # expand('{}/{{sample}}{{training_file}}'.format(ns_train_d), sample=config['samples'], training_file=ns_train_files),
        # expand('{}/{{gene}}/{{sample}}/{{data_file}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], data_file=gene_data),
        # expand('{}/{{gene}}/{{sample}}/reads_raw.{{extension}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], extension=['paf']),
        expand('{}/{{gene}}/{{sample}}/reads_raw.isoforms.{{extension}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], extension=['tsv']),
        # expand('{}/{{gene}}/{{sample}}/reads_raw.iterative_canonical_exons.{{extension}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], extension=['data', 'tsv', 'zeros_unaligned.tsv']),
        expand('{}/{{gene}}/{{sample}}/reads_raw.isoforms_plots.{{extension}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], extension=['pdf']),

rule freddie_make:
    input:
        'Makefile'
    output:
        config['exec']['freddie']
    conda:
        'freddie.env'
    threads:
        6
    shell:
        'make -j {{threads}} -C {}'.format(config['exec']['freddie'][0:config['exec']['freddie'].rfind('/')])

# rule nanosim_make:
#     output:
#         simulator = config['exec']['simulator'],
#         trainer   = config['exec']['read_analysis'],
#     params:
#         nanosim_d = '/'.join(config['exec']['simulator'].split('/')[:-2]),
#     conda:
#         'freddie.env'
#     shell:
#         'rm -r {{params.nanosim_d}}; git clone --branch {} {} {{params.nanosim_d}}'.format(config['nanosim']['tag'], config['nanosim']['git'])
#
rule reads_to_bam:
    input:
        reads  = lambda wildcards: config['samples'][wildcards.sample],
        target = lambda wildcards: config['references'][wildcards.target]
    output:
        bam = '{}/{{sample}}.{{target}}.sorted.bam'.format(mapped_d),
        bai = '{}/{{sample}}.{{target}}.sorted.bam.bai'.format(mapped_d),
        sam = '{}/{{sample}}.{{target}}.sorted.sam'.format(mapped_d),
    params:
        mapping_settings = lambda wildcards: config['mapping_settings'][wildcards.target]
    conda:
        'freddie.env'
    threads:
        32
    shell:
        'minimap2 -aY -x {params.mapping_settings} --eqx --MD -t {threads} {input.target} {input.reads} | '
        '  samtools sort -T {output.bam}.tmp -m 2G -@ {threads} -O bam - > {output.bam}; '
        'samtools index -b -@ {threads} {output.bam}; '
        'samtools view -h {output.bam} -o {output.sam}; '

rule get_gene_data:
    input:
        reads      = '{}/{{sample}}.dna.sorted.bam'.format(mapped_d),
        index      = '{}/{{sample}}.dna.sorted.bam.bai'.format(mapped_d),
        gtf        = config['annotations']['gtf'],
        genome     = config['references']['dna'],
        genome_fai = config['references']['dna_fai'],
        script     = config['exec']['gene_data'],
    params:
        out_dir=lambda wildcards: '{}/{}/{}'.format(genes_d, config['genes'][wildcards.gene], wildcards.sample)
    output:
        ['{}/{{gene}}/{{sample}}/{}'.format(genes_d, out_file) for out_file in gene_data]
    conda:
        'freddie.env'
    shell:
        '{input.script} -g {wildcards.gene} -t {input.gtf} -d {input.genome} -r {input.reads} -o {params.out_dir}'

rule process_polyA:
    input:
        reads  = '{}/{{gene}}/{{sample}}/{{read_type}}.fasta'.format(genes_d),
        script = config['exec']['polyA'],
    output:
        report  = '{}/{{gene}}/{{sample}}/{{read_type}}.polyA.tsv'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        'python {input.script} {input.reads} {output.report}'

rule trim_polyA:
    input:
        report  = '{}/{{gene}}/{{sample}}/{{read_type}}.polyA.tsv'.format(genes_d),
        script = config['exec']['trim'],
    output:
        reads  = '{}/{{gene}}/{{sample}}/{{read_type}}.polyA.trimmed.fasta'.format(genes_d),
        report  = '{}/{{gene}}/{{sample}}/{{read_type}}.polyA.trimmed.tsv'.format(genes_d),
    params:
        min_score = 20
    conda:
        'freddie.env'
    shell:
        'cat {input.report} | python {input.script} {output.report} {params.min_score} > {output.reads}'

rule freddie_align:
    input:
        reads  = '{}/{{gene}}/{{sample}}/{{read_type}}.polyA.trimmed.fasta'.format(genes_d),
        gene   = '{}/{{gene}}/{{sample}}/gene.fasta'.format(genes_d),
        script = config['exec']['freddie'],
    output:
        paf = '{}/{{gene}}/{{sample}}/{{read_type}}.paf'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} align -g {input.gene} -r {input.reads} > {output.paf}'

rule freddie_plot_simulated:
    input:
        paf             = '{}/{{gene}}/{{sample}}/reads_sim.oriented.paf'.format(genes_d),
        transcripts_tsv = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
        simulated_tsv   = '{}/{{gene}}/{{sample}}/reads_sim.oriented.tsv'.format(genes_d),
        script          = config['exec']['freddie'],
    output:
        dot             = '{}/{{gene}}/{{sample}}/reads_sim.oriented.dot'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} plot -p {input.paf} -a {input.transcripts_tsv} -s {input.simulated_tsv} > {output.dot}'

rule split_dot:
    input:
        dot           = '{}/{{gene}}/{{sample}}/reads_sim.oriented.dot'.format(genes_d),
        simulated_tsv = '{}/{{gene}}/{{sample}}/reads_sim.oriented.tsv'.format(genes_d),
        script        = config['exec']['split_dot'],
    output:
        done          = '{}/{{gene}}/{{sample}}/reads_sim.oriented.split_dot.done'.format(genes_d),
    params:
        prefix        = '{}/{{gene}}/{{sample}}/reads_sim.oriented.split.'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -d {input.dot} -t {input.simulated_tsv} -o {params.prefix}; '
        'touch {output.done} '

rule split_pdf:
    input:
        done   = '{}/{{gene}}/{{sample}}/reads_sim.oriented.split_dot.done'.format(genes_d),
    output:
        done   = '{}/{{gene}}/{{sample}}/reads_sim.oriented.split_pdf.done'.format(genes_d),
    params:
        prefix = '{}/{{gene}}/{{sample}}/reads_sim.oriented.split.'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        'for i in {params.prefix}*.dot; do prefix=${{i%.dot}}; cat $prefix.dot | dot -T pdf > $prefix.pdf; done; '
        'touch {output.done} '

rule dot_to_pdf:
    input:
        dot = '{}/{{gene}}/{{sample}}/reads_sim.oriented.dot'.format(genes_d),
    output:
        pdf = '{}/{{gene}}/{{sample}}/reads_sim.oriented.pdf'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        'cat {input.dot} | dot -T pdf > {output.pdf}'

rule disentangle:
    input:
        transcripts_tsv     = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
        paf                 = '{}/{{gene}}/{{sample}}/{{read_type}}.paf'.format(genes_d),
        canonical_exons_txt = '{}/{{gene}}/{{sample}}/{{read_type}}.canonical_exons.txt'.format(genes_d),
        script              = config['exec']['disentangle'],
    output:
        disentanglement = '{}/{{gene}}/{{sample}}/{{read_type}}.disentanglement.pdf'.format(genes_d),
        leaves = '{}/{{gene}}/{{sample}}/{{read_type}}.disentanglement.leaves.txt'.format(genes_d),
    params:
        out_prefix='{}/{{gene}}/{{sample}}/{{read_type}}.disentanglement'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -p {input.paf} -t {input.transcripts_tsv} -pk {input.canonical_exons_txt} -op {params.out_prefix}'

rule find_canonical_exon_iteratively:
    input:
        paf    = '{}/{{gene}}/{{sample}}/{{read_type}}.paf'.format(genes_d),
        script = config['exec']['find_canonical_exon_iteratively'],
    output:
        disentanglement = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.pdf'.format(genes_d),
        zeros_unaligned = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.zeros_unaligned.tsv'.format(genes_d),
        matrix          = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.data'.format(genes_d),
        exons           = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.tsv'.format(genes_d),
        read_names      = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.read_names.txt'.format(genes_d),
    params:
        out_prefix='{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -p {input.paf} -op {params.out_prefix}'

rule find_isoforms:
    input:
        matrix          = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.data'.format(genes_d),
        script          = config['exec']['find_isoforms'],
        exons           = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.tsv'.format(genes_d),
        unaligned_gaps  = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.tsv'.format(genes_d),
        zeros_unaligned = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.zeros_unaligned.tsv'.format(genes_d),
        read_names      = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.read_names.txt'.format(genes_d),
        polyA_report    = '{}/{{gene}}/{{sample}}/{{read_type}}.polyA.trimmed.tsv'.format(genes_d),
    output:
        isoforms        = '{}/{{gene}}/{{sample}}/{{read_type}}.isoforms.tsv'.format(genes_d),
    params:
        out_prefix      = '{}/{{gene}}/{{sample}}/{{read_type}}.isoforms'.format(genes_d),
        k               = 2,
        e               = 0.2,
        order_isoforms  = 'True',
        garbage_isoform = 'True',
        recycle_garbage = 'True',
        timeout         = 15,
        license         = config['gurobi_license'],
    threads:
        32
    conda:
        'freddie.env'
    shell:
        'export GRB_LICENSE_FILE={params.license}; '
        '{input.script} -d {input.matrix} -et {input.exons} '
        ' -names {input.read_names} -ptinfo {input.polyA_report} '
        ' -k {params.k} --garbage-isoform {params.garbage_isoform} --recycle-garbage {params.recycle_garbage} '
        ' -oi {params.order_isoforms} -e {params.e} -ug {input.zeros_unaligned}'
        ' -t {threads} -to {params.timeout} -op {params.out_prefix}'

rule plot_isoforms:
    input:
        paf      = '{}/{{gene}}/{{sample}}/{{read_type}}.paf'.format(genes_d),
        isoforms = '{}/{{gene}}/{{sample}}/{{read_type}}.isoforms.tsv'.format(genes_d),
        exons    = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.tsv'.format(genes_d),
        script   = config['exec']['plot_isoforms'],
    output:
        isoform_plot = '{}/{{gene}}/{{sample}}/{{read_type}}.isoforms_plots.pdf'.format(genes_d),
    params:
        out_prefix='{}/{{gene}}/{{sample}}/{{read_type}}.isoforms_plots'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -p {input.paf} -i {input.isoforms} -e {input.exons} -op {params.out_prefix}'
