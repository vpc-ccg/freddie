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
    'reads.filtered_out.fasta',
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
         # expand('{}/{{gene}}/{{sample}}/{{data_file}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], data_file=gene_data),
         expand('{}/{{gene}}/{{sample}}/reads.isoforms.{{extension}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], extension=['tsv']),
         expand('{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons.{{extension}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], extension=['data', 'tsv', 'zeros_unaligned.tsv']),
         expand('{}/{{gene}}/{{sample}}/reads.isoforms_plots.{{extension}}'.format(genes_d),   gene=config['genes'], sample=config['samples'], extension=['pdf']),

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

# rule cluster_paf_real:
#     input:
#         paf    = '{}/{{gene}}/{{sample}}/reads.paf'.format(genes_d),
#         script = config['exec']['clustering'],
#         transcripts_tsv = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
#     output:
#         cluster      = '{}/{{gene}}/{{sample}}/reads.cluster'.format(genes_d),
#         raw_coverage = '{}/{{gene}}/{{sample}}/reads.cluster.raw_coverage.txt'.format(genes_d),
#         log          = '{}/{{gene}}/{{sample}}/reads.cluster.log'.format(genes_d),
#         coverage_pdf = '{}/{{gene}}/{{sample}}/reads.cluster.coverage.pdf'.format(genes_d),
#     conda:
#         'freddie.env'
#     shell:
#         '{input.script} -p {input.paf} -o {output.cluster} -t {input.transcripts_tsv}'

rule find_canonical_exon:
    input:
        transcripts_tsv = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
        paf             = '{}/{{gene}}/{{sample}}/reads.paf'.format(genes_d),
        script          = config['exec']['find_canonical_exon'],
    output:
        canonical_exons_pdf = '{}/{{gene}}/{{sample}}/reads.canonical_exons.pdf'.format(genes_d),
        canonical_exons_txt = '{}/{{gene}}/{{sample}}/reads.canonical_exons.txt'.format(genes_d),
    params:
        out_prefix='{}/{{gene}}/{{sample}}/reads.canonical_exons'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -p {input.paf} -t {input.transcripts_tsv} -op {params.out_prefix}'

rule disentangle:
    input:
        transcripts_tsv     = '{}/{{gene}}/{{sample}}/transcripts.tsv'.format(genes_d),
        paf                 = '{}/{{gene}}/{{sample}}/reads.paf'.format(genes_d),
        canonical_exons_txt = '{}/{{gene}}/{{sample}}/reads.canonical_exons.txt'.format(genes_d),
        script              = config['exec']['disentangle'],
    output:
        disentanglement = '{}/{{gene}}/{{sample}}/reads.disentanglement.pdf'.format(genes_d),
        leaves = '{}/{{gene}}/{{sample}}/reads.disentanglement.leaves.txt'.format(genes_d),
    params:
        out_prefix='{}/{{gene}}/{{sample}}/reads.disentanglement'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -p {input.paf} -t {input.transcripts_tsv} -pk {input.canonical_exons_txt} -op {params.out_prefix}'

rule find_canonical_exon_iteratively:
    input:
        paf    = '{}/{{gene}}/{{sample}}/reads.paf'.format(genes_d),
        script = config['exec']['find_canonical_exon_iteratively'],
    output:
        # disentanglement = '{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons.pdf'.format(genes_d),
        zeros_unaligned = '{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons.zeros_unaligned.tsv'.format(genes_d),
        matrix          = '{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons.data'.format(genes_d),
        exons           = '{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons.tsv'.format(genes_d),
    params:
        out_prefix='{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -p {input.paf} -op {params.out_prefix}'

rule find_isoforms:
    input:
        matrix          = '{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons.data'.format(genes_d),
        script          = config['exec']['find_isoforms'],
        exons           = '{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons.tsv'.format(genes_d),
        unaligned_gaps  = '{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons.tsv'.format(genes_d),
        zeros_unaligned = '{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons.zeros_unaligned.tsv'.format(genes_d),
    output:
        log             = '{}/{{gene}}/{{sample}}/reads.isoforms.tsv'.format(genes_d),
    params:
        out_prefix      ='{}/{{gene}}/{{sample}}/reads.isoforms'.format(genes_d),
        k               = 20,
        e               = 0.2,
        order_isoforms  =True,
        timeout         =15,
    threads:
        32
    # conda:
    #     'freddie.env'
    shell:
        'PATH_OLD=$PATH; '
        'LD_LIBRARY_PATH_OLD=$LD_LIBRARY_PATH; '
        'export GUROBI_HOME="/home/borabi/docs/gurobi/gurobi811/linux64"; '
        'export PATH="${{GUROBI_HOME}}/bin:${{PATH}}"; '
        'export LD_LIBRARY_PATH="${{GUROBI_HOME}}/lib:${{LD_LIBRARY_PATH}}"; '
        'my_hostname=$(hostname); '
        'export GRB_LICENSE_FILE="/home/borabi/gurobi-licences/"$my_hostname".lic"; '
        'python2.7 {input.script} -k {params.k} -d {input.matrix} -op {params.out_prefix} -t {threads} -et {input.exons} -ug {input.zeros_unaligned} -e {params.e} -irp {params.order_isoforms} -to {params.timeout}; '
        'export PATH=$PATH_OLD; '
        'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_OLD; '

rule plot_isoforms:
    input:
        paf      = '{}/{{gene}}/{{sample}}/reads.paf'.format(genes_d),
        isoforms = '{}/{{gene}}/{{sample}}/reads.isoforms.tsv'.format(genes_d),
        exons    = '{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons.tsv'.format(genes_d),
        matrix   = '{}/{{gene}}/{{sample}}/reads.iterative_canonical_exons.data'.format(genes_d),
        script   = config['exec']['plot_isoforms'],
    output:
        isoform_plot = '{}/{{gene}}/{{sample}}/reads.isoforms_plots.pdf'.format(genes_d),
    params:
        out_prefix='{}/{{gene}}/{{sample}}/reads.isoforms_plots'.format(genes_d),
    conda:
        'freddie.env'
    shell:
        '{input.script} -p {input.paf} -i {input.isoforms} -e {input.exons} -d {input.matrix} -op {params.out_prefix}'
