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

genes_d    = '{}/genes'.format(outpath)
mapped_d   = '{}/mapped'.format(outpath)

gene_data=[
    'transcripts.tsv',
    'transcripts.fasta',
    'gene.fasta',
    'reads.fastq',
    'reads.paf',
]

rule all:
    input:
        expand('{}/{{sample}}.deSALT.sam'.format(mapped_d), sample=config['samples']),
        expand('{}/{{sample}}.deSALT.paf'.format(mapped_d), sample=config['samples']),
        expand('{}/{{gene}}/{{sample}}/{{data_file}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], data_file=gene_data),
        expand('{}/{{gene}}/{{sample}}/reads.segments.{{extension}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], extension=['txt', 'pdf']),
        # expand('{}/{{gene}}/{{sample}}/reads.isoforms.{{extension}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], extension=['tsv']),
        # expand('{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons.{{extension}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], extension=['data', 'tsv', 'zeros_unaligned.tsv']),
        # expand('{}/{{gene}}/{{sample}}/reads.isoforms_plots.{{extension}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], extension=['pdf']),

rule deSALT:
    input:
        index = config['references']['dna_desalt'],
        reads = lambda wildcards: config['samples'][wildcards.sample],
    output:
        sam = '{}/{{sample}}.deSALT.sam'.format(mapped_d),
    conda:
        'freddie.env'
    threads:
        42
    shell:
        'deSALT aln -t {threads} -x ont1d -o {output.sam} {input.index} {input.reads}'

rule sam_to_paf:
    input:
        script = config['exec']['sam_to_paf'],
        sam    = '{}/{{sample}}.deSALT.sam'.format(mapped_d),
    output:
        paf = '{}/{{sample}}.deSALT.paf'.format(mapped_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -s {input.sam} -p {output.paf}'

rule gene_data:
    input:
        reads  = lambda wildcards: config['samples'][wildcards.sample],
        gtf    = config['annotations']['gtf'],
        genome = config['references']['dna'],
        paf    = '{}/{{sample}}.deSALT.paf'.format(mapped_d),
        script = config['exec']['gene_data'],
    params:
        out_dir=lambda wildcards: '{}/{}/{}'.format(genes_d, config['genes'][wildcards.gene], wildcards.sample),
        padding=1000
    output:
        ['{}/{{gene}}/{{sample}}/{}'.format(genes_d, out_file) for out_file in gene_data]
    conda:
        'freddie.env'
    shell:
        '{input.script} --genome {input.genome} --gtf {input.gtf}'
        '  --fastqs {input.reads} --gene {wildcards.gene} --paf {input.paf}'
        '  --padding {params.padding} --output {params.out_dir}'

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

rule find_segments:
    input:
        transcripts_tsv = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
        paf             = '{}/{{gene}}/{{sample}}/{{read_type}}.paf'.format(genes_d),
        script          = config['exec']['segment'],
    output:
        pdf = '{}/{{gene}}/{{sample}}/{{read_type}}.segments.pdf'.format(genes_d),
        txt = '{}/{{gene}}/{{sample}}/{{read_type}}.segments.txt'.format(genes_d),
    threads:
        8
    params:
        out_prefix='{}/{{gene}}/{{sample}}/{{read_type}}.segments'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -p {input.paf} -t {input.transcripts_tsv} -op {params.out_prefix} -c {threads}'

rule find_isoforms:
    input:
        script          = config['exec']['find_isoforms'],
        matrix          = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.data'.format(genes_d),
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
        timeout         = config['gurobi']['timeout'],
        license         = config['gurobi']['license'],
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

rule score_clustering:
    input:
        script   = config['exec']['rand_index'],
        isoforms = '{}/{{gene}}/{{sample}}/{{read_type}}.isoforms.tsv'.format(genes_d),
        reads    = '{}/{{gene}}/{{sample}}/{{read_type}}.fasta'.format(genes_d),
    output:
        accuracy = '{}/{{gene}}/{{sample}}/{{read_type}}.rand'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        'python {input.script} -r {input.reads} -i {input.isoforms} -o {output.accuracy}'

rule plot_isoforms:
    input:
        paf             = '{}/{{gene}}/{{sample}}/{{read_type}}.paf'.format(genes_d),
        isoforms        = '{}/{{gene}}/{{sample}}/{{read_type}}.isoforms.tsv'.format(genes_d),
        exons           = '{}/{{gene}}/{{sample}}/{{read_type}}.iterative_canonical_exons.tsv'.format(genes_d),
        transcripts_tsv = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
        script          = config['exec']['plot_isoforms'],
    output:
        isoform_plot = '{}/{{gene}}/{{sample}}/{{read_type}}.isoforms_plots.pdf'.format(genes_d),
    params:
        out_prefix='{}/{{gene}}/{{sample}}/{{read_type}}.isoforms_plots'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -p {input.paf} -i {input.isoforms} -e {input.exons} -t {input.transcripts_tsv} -op {params.out_prefix}'
