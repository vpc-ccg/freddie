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

genes_d  = '{}/genes'.format(outpath)
mapped_d = '{}/mapped'.format(outpath)

gene_data=[
    'transcripts.tsv',
    'transcripts.fasta',
    'gene.fasta',
    'reads.fasta',
]
nanosim_read_analysis_files=[
    '_unaligned_length.pkl',
    '_aligned_region.pkl',
    '_aligned_reads.pkl',
    '_ht_length.pkl',
    '_match.hist',
    '_mis.hist',
    '_del.hist',
    '_ins.hist',
    '_first_match.hist',
    '_error_markov_model',
    '_ht_ratio.pkl',
    '.sam',
    '_match_markov_model',
    '_model_profile',
    '_error_rate.tsv',
]

rule all:
    input:
         # expand('{}/{{gene}}/{{sample}}/simulated_reads.oriented.split_pdf.done'.format(genes_d),   gene=config['genes'], sample=config['samples']),
         # expand('{}/{{gene}}/{{sample}}/simulated_reads.oriented.cluster'.format(genes_d),   gene=config['genes'], sample=config['samples']),
         expand('{}/{{gene}}/{{sample}}/reads.cluster.raw_coverage.meanshift_coverage.pdf'.format(genes_d),   gene=config['genes'], sample=config['samples']),
         # expand('{}/{{gene}}/{{sample}}/transcripts.disentanglement.txt'.format(genes_d),   gene=config['genes'], sample=config['samples']),

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

rule nanosim_make:
    output:
        config['exec']['read_analysis'],
        simulator = config['exec']['simulator'],
        zipped    = temp(config['url']['nanosim'].split('/')[-1]),
        unzipped  = temp(directory('NanoSim-{}/'.format(config['url']['nanosim'].split('/')[-1].split('.')[0]))),
    conda:
        'freddie.env'
    shell:
        '''
        wget {};
        unzip {{output.zipped}};
        cp {{output.unzipped}}/* extern/nanosim/ -r;
        patch {{output.simulator}} extern/nanosim.patch
        '''.format(config['url']['nanosim'])

rule minimap2_map:
    input:
        genome=config['references']['genome'],
        reads=lambda wildcards: config['samples'][wildcards.sample],
    output:
        temp('{}/{{sample}}.sam'.format(mapped_d))
    conda:
        'freddie.env'
    threads:
        32
    shell:
        'minimap2 -aY -x splice --eqx -t {threads} {input.genome} {input.reads} > {output}'

rule samtools_sort:
    input:
        '{}/{{sample}}.sam'.format(mapped_d),
    output:
        bam='{}/{{sample}}.sorted.bam'.format(mapped_d),
    conda:
        'freddie.env'
    threads:
        32
    shell:
        'samtools sort -T {}/{{wildcards.sample}}.tmp -m 2G -@ {{threads}} -O bam {{input}} > {{output.bam}}'.format(mapped_d)

rule samtools_index:
    input:
        bam='{}/{{sample}}.sorted.bam'.format(mapped_d),
    output:
        index='{}/{{sample}}.sorted.bam.bai'.format(mapped_d),
    conda:
        'freddie.env'
    threads:
        32
    shell:
        'samtools index -b -@ {threads} {input.bam}'

rule get_gene_data:
    input:
        reads  =  '{}/{{sample}}.sorted.bam'.format(mapped_d),
        index  =  '{}/{{sample}}.sorted.bam.bai'.format(mapped_d),
        gtf    =  config['annotations']['gtf'],
        genome =  config['references']['genome'],
        script =  config['exec']['gene_data'],
    params:
        out_dir=lambda wildcards: '{}/{}/{}'.format(genes_d, config['genes'][wildcards.gene], wildcards.sample)
    output:
        ['{}/{{gene}}/{{sample}}/{}'.format(genes_d, out_file) for out_file in gene_data]
    conda:
        'freddie.env'
    shell:
        '{input.script} -g {wildcards.gene} -t {input.gtf} -d {input.genome} -r {input.reads} -o {params.out_dir}'

rule freddie_align_real:
    input:
        reads='{}/{{gene}}/{{sample}}/reads.fasta'.format(genes_d),
        gene = '{}/{{gene}}/{{sample}}/gene.fasta'.format(genes_d),
        script = config['exec']['freddie'],
    output:
        paf = '{}/{{gene}}/{{sample}}/reads.paf'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} align -g {input.gene} -r {input.reads} > {output.paf}'

rule train_nanosim:
    input:
        transcripts='{}/{{gene}}/{{sample}}/transcripts.fasta'.format(genes_d),
        reads='{}/{{gene}}/{{sample}}/reads.fasta'.format(genes_d),
        script =  config['exec']['read_analysis'],
    params:
        out_prefix='{}/{{gene}}/{{sample}}/nanosim/training'.format(genes_d),
    output:
        ['{}/{{gene}}/{{sample}}/nanosim/training{}'.format(genes_d, training_file) for  training_file in nanosim_read_analysis_files]
    conda:
        'freddie.env'
    threads:
        32
    shell:
        '{input.script} -i {input.reads} -r {input.transcripts} -t {threads} -o {params.out_prefix}'

rule run_nanosim:
    input:
        ['{}/{{gene}}/{{sample}}/nanosim/training{}'.format(genes_d, training_file) for  training_file in nanosim_read_analysis_files],
        transcript_tsv   = '{}/{{gene}}/{{sample}}/{}'.format(genes_d, 'transcripts.tsv'),
        transcript_fasta = '{}/{{gene}}/{{sample}}/transcripts.fasta'.format(genes_d),
        nanosim          = config['exec']['simulator'],
        script           = config['exec']['run_nanosim']
    params:
        train_prefix='{}/{{gene}}/{{sample}}/nanosim/training'.format(genes_d),
        intermediate_directory='{}/{{gene}}/{{sample}}/nanosim/'.format(genes_d),
        read_count=100,
        distribution='normal'
    output:
        simulated_tsv='{}/{{gene}}/{{sample}}/simulated_reads.oriented.tsv'.format(genes_d),
        oriented_reads = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.fasta'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -tt {input.transcript_tsv} -tf {input.transcript_fasta}'
        ' -ns {input.nanosim} -tr {params.train_prefix}'
        ' -d {params.intermediate_directory} -c {params.read_count} -f {params.distribution}'
        ' -or {output.oriented_reads} -ot {output.simulated_tsv}'

rule freddie_align:
    input:
        reads = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.fasta'.format(genes_d),
        gene = '{}/{{gene}}/{{sample}}/gene.fasta'.format(genes_d),
        script = config['exec']['freddie'],
    output:
        paf = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.paf'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} align -g {input.gene} -r {input.reads} > {output.paf}'

rule freddie_plot:
    input:
        paf = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.paf'.format(genes_d),
        transcripts_tsv = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
        simulated_tsv='{}/{{gene}}/{{sample}}/simulated_reads.oriented.tsv'.format(genes_d),
        script = config['exec']['freddie'],
    output:
        dot = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.dot'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} plot -p {input.paf} -a {input.transcripts_tsv} -s {input.simulated_tsv} > {output.dot}'

rule split_dot:
    input:
        dot = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.dot'.format(genes_d),
        simulated_tsv='{}/{{gene}}/{{sample}}/simulated_reads.oriented.tsv'.format(genes_d),
        script = config['exec']['split_dot'],
    output:
        done = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.split_dot.done'.format(genes_d),
    params:
        prefix = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.split.'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -d {input.dot} -t {input.simulated_tsv} -o {params.prefix}; '
        'touch {output.done} '

rule split_pdf:
    input:
        done = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.split_dot.done'.format(genes_d),
    output:
        done = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.split_pdf.done'.format(genes_d),
    params:
        prefix = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.split.'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        'for i in {params.prefix}*.dot; do prefix=${{i%.dot}}; cat $prefix.dot | dot -T pdf > $prefix.pdf; done; '
        'touch {output.done} '

rule dot_to_pdf:
    input:
        dot = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.dot'.format(genes_d),
    output:
        pdf = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.pdf'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        'cat {input.dot} | dot -T pdf > {output.pdf}'

rule cluster_paf_real:
    input:
        paf    = '{}/{{gene}}/{{sample}}/reads.paf'.format(genes_d),
        script = config['exec']['clustering'],
        transcripts_tsv = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
    output:
        cluster      = '{}/{{gene}}/{{sample}}/reads.cluster'.format(genes_d),
        raw_coverage = '{}/{{gene}}/{{sample}}/reads.cluster.raw_coverage.txt'.format(genes_d),
        log          = '{}/{{gene}}/{{sample}}/reads.cluster.log'.format(genes_d),
        coverage_pdf = '{}/{{gene}}/{{sample}}/reads.cluster.coverage.pdf'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -p {input.paf} -o {output.cluster} -t {input.transcripts_tsv}'

rule meanshift_real:
    input:
        transcripts_tsv = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
        paf          = '{}/{{gene}}/{{sample}}/reads.paf'.format(genes_d),
        raw_coverage = '{}/{{gene}}/{{sample}}/reads.cluster.raw_coverage.txt'.format(genes_d),
        script       = config['exec']['meanshift'],
    output:
        meanshift_coverage_pdf = '{}/{{gene}}/{{sample}}/reads.cluster.raw_coverage.meanshift_coverage.pdf'.format(genes_d),
    params:
        bin_size=0.025,
        window_size=1,
    conda:
        'freddie.env'
    shell:
        '{input.script} -c {input.raw_coverage} -bs {params.bin_size} -w {params.window_size} -p {input.paf} -t {input.transcripts_tsv}'

rule disentangle:
    input:
        transcripts_tsv = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
        script          = config['exec']['disentangle'],
    output:
        disentanglement = '{}/{{gene}}/{{sample}}/transcripts.disentanglement.txt'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -t {input.transcripts_tsv} -o {output.disentanglement}'


# rule cluster_paf:
#     input:
#         paf    = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.paf'.format(genes_d),
#         script = config['exec']['clustering'],
#         transcripts_tsv = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
#     output:
#         cluster      = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.cluster'.format(genes_d),
#         raw_coverage = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.cluster.raw_coverage.txt'.format(genes_d),
#         log          = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.cluster.log'.format(genes_d),
#         coverage_pdf = '{}/{{gene}}/{{sample}}/simulated_reads.oriented.cluster.coverage.pdf'.format(genes_d),
#     conda:
#         'freddie.env'
#     shell:
#         '{input.script} -p {input.paf} -o {output.cluster} -t {input.transcripts_tsv}'
