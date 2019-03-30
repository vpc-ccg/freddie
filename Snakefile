configfile: 'config.yaml'

outpath = config['outpath'].rstrip('/')
if outpath in ['.', './'] :
    outpath = ''

genes_d  = '{}genes'.format(outpath)
mapped_d = '{}mapped'.format(outpath)

rule all:
    input:
         expand(directory(genes_d+'/{gene}/{sample}'), gene=config['genes'], sample=config['samples']),
         # expand(genes_d+'/{gene}/{sample}/{gene_file}', gene=config['genes'], sample=config['samples'], gene_file=['transcripts.tsv', 'transcripts.fasta', 'gene.fasta', 'reads.fasta']),

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
        gene=lambda wildcards: config['genes'][wildcards.gene]
    output:
        directory=directory(genes_d+'/{gene}/{sample}'),
    conda:
        'freddie.env'
    shell:
        '{input.script} -g {params.gene} -t {input.gtf} -d {input.genome} -r {input.reads} -o {output.directory}'











print()
