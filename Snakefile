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
        # expand('{}/{{sample}}.deSALT.sam'.format(mapped_d), sample=config['samples']),
        # expand('{}/{{sample}}.deSALT.paf'.format(mapped_d), sample=config['samples']),
        # expand('{}/{{gene}}/{{sample}}/{{data_file}}'.format(genes_d),
        #     gene=config['genes'], sample=config['samples'], data_file=gene_data),
        expand('{}/{{gene}}/{{sample}}/reads.segments.{{extension}}'.format(genes_d),
            gene=config['genes'], sample=config['samples'], extension=['txt', 'pdf', 'data', 'gaps', 'names']),
        expand('{}/{{gene}}/{{sample}}/reads.isoforms.{{recycle_model}}.{{extension}}'.format(genes_d),
            gene=config['genes'], sample=config['samples'], recycle_model=config['gurobi']['recycle_models'], extension=['tsv']),
        expand('{}/{{gene}}/{{sample}}/reads.isoforms.{{recycle_model}}.plots.pdf'.format(genes_d),
            gene=config['genes'], sample=config['samples'], recycle_model=config['gurobi']['recycle_models'],),

rule deSALT:
    input:
        index = config['references']['dna_desalt'],
        reads = lambda wildcards: config['samples'][wildcards.sample],
    output:
        sam = protected('{}/{{sample}}.deSALT.sam'.format(mapped_d)),
    params:
        temp = '{}/{{sample}}.deSALT.temp.'.format(mapped_d),
    conda:
        'freddie.env'
    threads:
        20
    shell:
        'deSALT aln -t {threads} -x ont1d --temp-file-perfix {params.temp} -o {output.sam} {input.index} {input.reads}'

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

# rule process_polyA:
#     input:
#         reads  = '{}/{{gene}}/{{sample}}/reads.fasta'.format(genes_d),
#         script = config['exec']['polyA'],
#     output:
#         report  = '{}/{{gene}}/{{sample}}/reads.polyA.tsv'.format(genes_d),
#     conda:
#         'freddie.env'
#     shell:
#         'python {input.script} {input.reads} {output.report}'
#
# rule trim_polyA:
#     input:
#         report  = '{}/{{gene}}/{{sample}}/reads.polyA.tsv'.format(genes_d),
#         script = config['exec']['trim'],
#     output:
#         reads  = '{}/{{gene}}/{{sample}}/reads.polyA.trimmed.fasta'.format(genes_d),
#         report  = '{}/{{gene}}/{{sample}}/reads.polyA.trimmed.tsv'.format(genes_d),
#     params:
#         min_score = 20
#     conda:
#         'freddie.env'
#     shell:
#         'cat {input.report} | python {input.script} {output.report} {params.min_score} > {output.reads}'

rule find_segments:
    input:
        transcripts_tsv = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
        paf             = '{}/{{gene}}/{{sample}}/reads.paf'.format(genes_d),
        script          = config['exec']['segment'],
    output:
        [protected('{}/{{gene}}/{{sample}}/reads.segments.{}'.format(genes_d,ext)) for ext in ['pdf','txt','data','gaps','names']]
    threads:
        8
    params:
        out_prefix='{}/{{gene}}/{{sample}}/reads.segments'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -p {input.paf} -t {input.transcripts_tsv} -op {params.out_prefix} -c {threads}'

rule find_isoforms:
    input:
        script  = config['exec']['find_isoforms'],
        segs    = '{}/{{gene}}/{{sample}}/reads.segments.txt'.format(genes_d),
        data    = '{}/{{gene}}/{{sample}}/reads.segments.data'.format(genes_d),
        gaps    = '{}/{{gene}}/{{sample}}/reads.segments.gaps'.format(genes_d),
        names   = '{}/{{gene}}/{{sample}}/reads.segments.names'.format(genes_d),
    output:
        isoforms = protected('{}/{{gene}}/{{sample}}/reads.isoforms.{{recycle_model}}.tsv'.format(genes_d)),
    params:
        out_prefix = '{}/{{gene}}/{{sample}}/reads.isoforms.{{recycle_model}}'.format(genes_d),
        epsilon    = 0.2,
        timeout    = config['gurobi']['timeout'],
        license    = config['gurobi']['license'],
    threads:
        32
    conda:
        'freddie.env'
    shell:
        'export GRB_LICENSE_FILE={params.license}; '
        '{input.script}'
        ' -rm    {wildcards.recycle_model} '
        ' -st    {input.segs} '
        ' -sd    {input.data} '
        ' -ug    {input.gaps}'
        ' -e     {params.epsilon} '
        ' -names {input.names} '
        ' -t     {threads}'
        ' -to    {params.timeout} '
        ' -op    {params.out_prefix}'

rule plot_isoforms:
    input:
        script      = config['exec']['plot_isoforms'],
        transcripts = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
        segs        = '{}/{{gene}}/{{sample}}/reads.segments.txt'.format(genes_d),
        isoforms    = '{}/{{gene}}/{{sample}}/reads.isoforms.{{recycle_model}}.tsv'.format(genes_d),
    output:
        isoform_pdf = protected('{}/{{gene}}/{{sample}}/reads.isoforms.{{recycle_model}}.plots.pdf'.format(genes_d)),
    params:
        out_prefix  = '{}/{{gene}}/{{sample}}/reads.isoforms.{{recycle_model}}.plots'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -t {input.transcripts} -s {input.segs} -i {input.isoforms} -op {params.out_prefix}'
