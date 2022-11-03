configfile: 'config.yaml'

for k in config.keys():
    if not k.startswith('override_'):
        continue
    keys = k[len('override_'):].split('_')
    top_dict = eval('config{}'.format(''.join(['["{}"]'.format(x) for x in keys[:-1]])))
    assert keys[-1] in top_dict
    top_dict[keys[-1]]=config[k]

outpath = config['outpath'].rstrip('/')

output_d = '{}/results'.format(outpath)
logs_d   = '{}/logs'.format(outpath)
benchmark_d   = '{}/benchmark'.format(outpath)

rule all:
    input:
        expand('{}/{{sample}}/{{sample}}.sorted.bam'.format(output_d), sample=config['samples']),
        expand('{}/{{sample}}/{{run_mode}}/freddie.split'.format(output_d),         sample=config['samples'], run_mode=config['run_modes']),
        expand('{}/{{sample}}/{{run_mode}}/freddie.segment'.format(output_d),       sample=config['samples'], run_mode=config['run_modes']),
        expand('{}/{{sample}}/{{run_mode}}/freddie.cluster'.format(output_d),       sample=config['samples'], run_mode=config['run_modes']),
        expand('{}/{{sample}}/{{run_mode}}/freddie.isoforms.gtf'.format(output_d),  sample=config['samples'], run_mode=config['run_modes']),

rule minimap2:
    input:
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        genome = lambda wildcards: config['references'][config['samples'][wildcards.sample]['ref']]['genome'],
    output:
        bam=protected('{}/{{sample}}/{{sample}}.sorted.bam'.format(output_d)),
        bai=protected('{}/{{sample}}/{{sample}}.sorted.bam.bai'.format(output_d)),
    conda:
        'envs/minimap2.yml'
    threads:
        64
    resources:
        mem  = "128G",
        time = 1439,
    shell:
        'minimap2 -a -x splice -t {threads} {input.genome} {input.reads} | '
        '  samtools sort -T {output.bam}.tmp -m 2G -@ {threads} -O bam - > {output.bam} && '
        '  samtools index {output.bam} '

rule split:
    input:
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        bam = '{}/{{sample}}/{{sample}}.sorted.bam'.format(output_d),
    output:
        split = directory('{}/{{sample}}/{{run_mode}}/freddie.split'.format(output_d)),
    params:
        script = config['exec']['split'],
        run_mode = lambda wildcards: config['run_modes'][wildcards.run_mode]['split']
    conda:
        'envs/freddie.yml'
    threads:
        32
    resources:
        mem  = "16G",
        time = 359,
    benchmark:
        f'{benchmark_d}/split/{{sample}}.{{run_mode}}.tsv'
    shell:
        '{params.script} -b {input.bam} -r {input.reads} -o {output.split} -t {threads} {params.run_mode}'

rule segment:
    input:
        split  = '{}/{{sample}}/{{run_mode}}/freddie.split'.format(output_d),
    output:
        segment = directory('{}/{{sample}}/{{run_mode}}/freddie.segment'.format(output_d)),
    params:
        script = config['exec']['segment'],
        run_mode = lambda wildcards: config['run_modes'][wildcards.run_mode]['segment']
    conda:
        'envs/freddie.yml'
    threads:
        32
    resources:
        mem  = "32G",
        time = 599,
    benchmark:
        f'{benchmark_d}/segment/{{sample}}.{{run_mode}}.tsv'
    shell:
        '{params.script} -s {input.split} -o {output.segment} -t {threads} {params.run_mode}'

rule cluster:
    input:
        license = config['gurobi']['license'],
        segment = '{}/{{sample}}/{{run_mode}}/freddie.segment'.format(output_d),
    output:
        cluster = directory('{}/{{sample}}/{{run_mode}}/freddie.cluster'.format(output_d)),
        logs    = directory('{}/{{sample}}/{{run_mode}}/freddie.cluster_logs'.format(output_d)),
        log     = '{}/{{sample}}/{{run_mode}}/freddie.cluster.log'.format(logs_d),
    conda:
        'envs/freddie.yml'
    params:
        script  = config['exec']['cluster'],
        timeout = config['gurobi']['timeout'],
        run_mode = lambda wildcards: config['run_modes'][wildcards.run_mode]['cluster']
    conda:
        'envs/freddie.yml'
    threads:
        32
    resources:
        mem  = "64G",
        time = 999,
    benchmark:
        f'{benchmark_d}/cluster/{{sample}}.{{run_mode}}.tsv'
    shell:
        'export GRB_LICENSE_FILE={input.license}; '
        '{params.script} -s {input.segment} -o {output.cluster} -l {output.logs} -t {threads} -to {params.timeout} {params.run_mode} > {output.log}'

rule isoforms:
    input:
        script  = config['exec']['isoforms'],
        split   = '{}/{{sample}}/freddie.split'.format(output_d),
        cluster = '{}/{{sample}}/freddie.cluster'.format(output_d),
    output:
        isoforms = protected('{}/{{sample}}/freddie.isoforms.gtf'.format(output_d)),
    conda:
        'envs/freddie.yml'
    threads:
        8
    resources:
        mem  = "16G",
        time = 359,
    shell:
        '{input.script} -s {input.split} -c {input.cluster} -o {output.isoforms} -t {threads}'