configfile: 'config.yaml'

outpath = config['outpath'].rstrip('/')
if outpath in ['.', './'] :
    outpath = ''

genes_d  = '{}genes'.format(outpath)
mapped_d = '{}mapped'.format(outpath)
gene_data=[
    'transcripts.tsv',
    'transcripts.fasta',
    'gene.fasta',
    'reads.fasta',
]
nanosim_read_analysis_files=[
    '_aligned_region.pkl',
    '_aligned_reads.pkl',
    '_ht_length.pkl',
    '_besthit.maf/sam',
    '_match.hist/training_mis.hist/training_del.hist/training_ins.hist',
    '_first_match.hist',
    '_error_markov_model',
    '_ht_ratio.pkl',
    '_maf/sam',
    '_match_markov_model',
    '_model_profile',
    '_processed.maf',
    '_unaligned_length.pkl',
    '_error_rate.tsv',
]
nanosim_simulator_files=[
    '_reads.fasta',
    '.log',
]

rule all:
    input:
         expand(genes_d+'/{gene}/{sample}/{out_file}', gene=config['genes'], sample=config['samples'], out_file=gene_data),
         expand(genes_d+'/{gene}/{sample}/simulated{out_file}', gene=config['genes'], sample=config['samples'], out_file=nanosim_simulator_files),
         expand(genes_d+'/{gene}/{sample}/simulated_reads.tsv', gene=config['genes'], sample=config['samples']),
         expand(genes_d+'/{gene}/{sample}/simulated_reads.pdf', gene=config['genes'], sample=config['samples']),

rule freddie_make:
    input:
        'Makefile'
    output:
        config['exec']['freddie']
    shell:
        'make -Bj'

rule nanosim_make:
    output:
        config['exec']['read_analysis'],
        config['exec']['simulator'],
        zipped=temp(config['url']['nanosim'].split('/')[-1]),
        unzipped=temp(directory('NanoSim-{}/'.format(config['url']['nanosim'].split('/')[-1].split('.')[0]))),
    conda:
        'freddie.env'
    shell:
        'wget ' + config['url']['nanosim'] + '; '
        'unzip {output.zipped}; '
        'cp {output.unzipped}/* extern/nanosim/ -r '

rule minimap2_map:
    input:
        genome=config['references']['genome'],
        reads=lambda wildcards: config['samples'][wildcards.sample],
    output:
        temp(mapped_d+'/{sample}.sam')
    conda:
        'freddie.env'
    threads:
        32
    shell:
        'minimap2 -aY -x splice -t -eqx -x -t {threads} {input.genome} {input.reads} > {output}'

rule samtools_sort:
    input:
        mapped_d+'/{sample}.sam'
    output:
        bam=mapped_d+'/{sample}.sorted.bam',
        index=mapped_d+'/{sample}.sorted.bam.bai'
    conda:
        'freddie.env'
    threads: 32
    shell:
        'samtools sort -T '+ mapped_d +'/{wildcards.sample}.tmp -m 2G -@ {threads} -O bam {input} > {output.bam} ; '
        'samtools index -m 2G -@ {threads} {output.bam} ; '

rule get_gene_data:
    input:
        reads  =  mapped_d+'/{sample}.sorted.bam',
        gtf    =  config['annotations']['gtf'],
        genome =  config['references']['genome'],
        script =  config['exec']['gene_data'],
    params:
        out_dir=lambda wildcards: '{}/{}/{}'.format(genes_d, config['genes'][wildcards.gene], config['samples'][wildcards.sample])
    output:
        [genes_d+'/{gene}/{sample}/'+out_file for out_file in gene_data]
    conda:
        'freddie.env'
    shell:
        '{input.script} -g {wildcards.gene} -t {input.gtf} -d {input.genome} -r {input.reads} -o {params.out_dir}'

rule nanosim_read_analysis:
    input:
        transcripts=genes_d+'/{gene}/{sample}/transcripts.fasta',
        reads=genes_d+'/{gene}/{sample}/reads.fasta',
        script =  config['exec']['read_analysis'],
    params:
        out_prefix=genes_d+'/{gene}/{sample}/training',
    output:
        [genes_d+'/{gene}/{sample}/training'+training_file for  training_file in nanosim_read_analysis_files]
    conda:
        'freddie.env'
    threads:
        32
    shell:
        '{input.script} -i {input.reads} -r {input.transcripts} -r {input.reads} -t {threads} -o {params.out_prefix}'

rule nanosim_simulate:
    input:
        [genes_d+'/{gene}/{sample}/training'+training_file for  training_file in nanosim_read_analysis_files],
        transcripts=genes_d+'/{gene}/{sample}/transcripts.fasta',
        script =  config['exec']['simulator'],
    params:
        in_prefix=genes_d+'/{gene}/{sample}/training',
        out_prefix=genes_d+'/{gene}/{sample}/simulated',
        read_count=10
    output:
        [genes_d+'/{gene}/{sample}/simulated'+simulation_file for  simulation_file in nanosim_simulator_files]
    conda:
        'freddie.env'
    shell:
        '{input.script} linear -r {input.transcripts} -c {params.in_prefix} -o {params.out_prefix} -n {params.read_count}'

rule get_nanosim_tsv:
    input:
        reads           = genes_d+'/{gene}/{sample}/simulated_reads.fasta',
        transcripts_tsv = genes_d+'/{gene}/{sample}/transcripts.tsv',
        script          = config['exec']['nanosim_tsv'],
    output:
        simulated_tsv=genes_d+'/{gene}/{sample}/simulated_reads.tsv',
    conda:
        'freddie.env'
    shell:
        '{input.script} -nsr {input.reads} -t {input.transcripts_tsv} -o {output.simulated_tsv}'

rule freddie_align:
    input:
        reads = genes_d+'/{gene}/{sample}/simulated_reads.fasta',
        gene = genes_d+'/{gene}/{sample}/gene.fasta',
        script = config['exec']['freddie'],
    output:
        paf = genes_d+'/{gene}/{sample}/simulated_reads.paf',
    conda:
        'freddie.env'
    shell:
        '{input.script} align -g {input.gene} -r {input.reads} > {output.paf}'

rule freddie_plot:
    input:
        paf = genes_d+'/{gene}/{sample}/simulated_reads.paf',
        transcripts_tsv = genes_d+'/{gene}/{sample}/transcripts.tsv',
        simulated_tsv = genes_d+'/{gene}/{sample}/simulated_reads.tsv',
        script = config['exec']['freddie'],
    output:
        dot = genes_d+'/{gene}/{sample}/simulated_reads.dot',
    conda:
        'freddie.env'
    shell:
        '{input.script} plot -p {input.paf} -a {input.transcripts_tsv} -s {input.simulated_tsv} > {output.dot}'

rule dot_to_pdf:
    input:
        dot = genes_d+'/{gene}/{sample}/simulated_reads.dot',
    output:
        pdf = genes_d+'/{gene}/{sample}/simulated_reads.pdf',
    conda:
        'freddie.env'
    shell:
        'cat {input.dot} | dot -T pdf > {output.pdf}'
