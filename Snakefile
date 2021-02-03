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
mapped_d = '{}/mapped'.format(outpath)
logs_d   = '{}/logs'.format(outpath)

rule all:
    input:
        expand('{}/{{sample}}/{{sample}}.bam'.format(mapped_d),       sample=config['samples']),
        expand('{}/{{sample}}/freddie.split'.format(output_d),     sample=config['samples']),
        expand('{}/{{sample}}/freddie.segment'.format(output_d),   sample=config['samples']),
        expand('{}/{{sample}}/freddie.cluster'.format(output_d),   sample=config['samples']),
        expand('{}/{{sample}}/freddie.isoforms.gtf'.format(output_d),  sample=config['samples']),

rule align:
    input:
        script = config['exec']['align'],
        index  = config['references']['dna_desalt'],
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
    output:
        sam=protected('{}/{{sample}}/{{sample}}.sam'.format(mapped_d)),
        bam=protected('{}/{{sample}}/{{sample}}.bam'.format(mapped_d)),
        bai=protected('{}/{{sample}}/{{sample}}.bam.bai'.format(mapped_d)),
    params:
        seq_type = lambda wildcards: config['samples'][wildcards.sample]['seq_type'],
    conda:
        'environment.env'
    threads:
        32
    shell:
        'py/freddie_align.py -r {input.reads} -i {input.index}/ -s {params.seq_type} -o {output.sam} -t {threads} &&'
        '  samtools sort -T {output.bam}.tmp -m 2G -@ {threads} -O bam {output.sam} > {output.bam} && '
        '  samtools index {output.bam} '

rule split:
    input:
        script = config['exec']['split'],
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        bam = '{}/{{sample}}/{{sample}}.bam'.format(mapped_d),
    output:
        split = directory('{}/{{sample}}/freddie.split'.format(output_d)),
    conda:
        'environment.env'
    threads:
        8
    shell:
        '{input.script} -b {input.bam} -r {input.reads} -o {output.split} -t {threads}'

rule segment:
    input:
        script = config['exec']['segment'],
        split  = '{}/{{sample}}/freddie.split'.format(output_d),
    output:
        segment = directory('{}/{{sample}}/freddie.segment'.format(output_d)),
    conda:
        'environment.env'
    threads:
        32
    shell:
        '{input.script} -s {input.split} -o {output.segment} -t {threads}'

rule cluster:
    input:
        license = config['gurobi']['license'],
        script  = config['exec']['cluster'],
        segment = '{}/{{sample}}/freddie.segment'.format(output_d),
    output:
        cluster = directory('{}/{{sample}}/freddie.cluster'.format(output_d)),
        logs    = directory('{}/{{sample}}/freddie.cluster_logs'.format(output_d)),
        log     = '{}/{{sample}}/freddie.cluster.log'.format(logs_d),
    conda:
        'environment.env'
    params:
        timeout = config['gurobi']['timeout'],
    threads:
        32
    shell:
        'export GRB_LICENSE_FILE={input.license}; '
        '{input.script} -s {input.segment} -o {output.cluster} -l {output.logs} -t {threads} -to {params.timeout} > {output.log}'

rule isoforms:
    input:
        script  = config['exec']['isoforms'],
        segment = '{}/{{sample}}/freddie.segment'.format(output_d),
        cluster = '{}/{{sample}}/freddie.cluster'.format(output_d),
    output:
        isoforms = protected('{}/{{sample}}/freddie.isoforms.gtf'.format(output_d)),
    shell:
        '{input.script} -s {input.segment} -c {input.cluster} -o {output.isoforms}'