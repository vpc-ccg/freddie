configfile: 'config.yaml'

for k in config.keys():
    if not k.startswith('override_'):
        continue
    keys = k[len('override_'):].split('_')
    top_dict = eval('config{}'.format(''.join(['["{}"]'.format(x) for x in keys[:-1]])))
    assert keys[-1] in top_dict
    top_dict[keys[-1]]=config[k]

def make_slurm():
    import os
    os.makedirs('slurm'.format(outpath), mode=0o777, exist_ok=True)

outpath = config['outpath'].rstrip('/')
make_slurm()

output_d = '{}/results'.format(outpath)
plots_d  = '{}/plots'.format(outpath)
mapped_d = '{}/mapped'.format(outpath)
logs_d   = '{}/logs'.format(outpath)

rule all:
    input:
        # expand('{}/{{sample}}.ref_cov.tsv'.format(output_d),           sample=config['samples']),
        # expand('{}/{{sample}}.ref_beds'.format(output_d),             sample=config['samples']),
        expand('{}/freddie.{{sample}}.sam'.format(mapped_d),           sample=config['samples']),
        expand('{}/freddie.{{sample}}.split.tsv'.format(output_d),     sample=config['samples']),
        expand('{}/freddie.{{sample}}.segment.tsv'.format(output_d),   sample=config['samples']),
        expand('{}/freddie.{{sample}}.cluster.tsv'.format(output_d),   sample=config['samples']),
        expand('{}/freddie.{{sample}}.isoforms.gtf'.format(output_d),  sample=config['samples']),
        expand('{}/freddie.{{sample}}.beds'.format(output_d),         sample=config['samples']),
        expand('{}/freddie.{{sample}}.intersect.tsv'.format(output_d), sample=config['samples']),
        expand('{}/freddie.{{sample}}.seqpare.tsv'.format(output_d),   sample=config['samples']),
        expand('{}/freddie.{{sample}}.pairings.tsv'.format(output_d),  sample=config['samples']),
        # expand('{}/freddie.{{sample}}.plot'.format(plots_d),         sample=config['samples']),

rule align:
    input:
        script = config['exec']['align'],
        index  = config['references']['dna_desalt'],
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
    output:
        sam = protected('{}/freddie.{{sample}}.sam'.format(mapped_d)),
    params:
        seq_type = lambda wildcards: config['samples'][wildcards.sample]['seq_type'],
    conda:
        'environment.env'
    threads:
        32
    shell:
        '{input.script} -r {input.reads} -i {input.index}/ -s {params.seq_type} -o {output.sam} -t {threads}'

rule split:
    input:
        script = config['exec']['split'],
        sam    = '{}/freddie.{{sample}}.sam'.format(mapped_d),
    output:
        split = protected('{}/freddie.{{sample}}.split.tsv'.format(output_d)),
    conda:
        'environment.env'
    shell:
        '{input.script} -s {input.sam} -o {output.split}'

rule segment:
    input:
        script = config['exec']['segment'],
        split  = '{}/freddie.{{sample}}.split.tsv'.format(output_d),
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
    output:
        segment = protected('{}/freddie.{{sample}}.segment.tsv'.format(output_d)),
    conda:
        'environment.env'
    threads:
        32
    shell:
        '{input.script} -s {input.split} -r {input.reads} -o {output.segment} -t {threads}'

rule cluster:
    input:
        license = config['gurobi']['license'],
        script  = config['exec']['cluster'],
        segment ='{}/freddie.{{sample}}.segment.tsv'.format(output_d),
    output:
        cluster = protected('{}/freddie.{{sample}}.cluster.tsv'.format(output_d)),
        logs    = directory('{}/freddie.{{sample}}.cluster'.format(logs_d)),
        log     = '{}/freddie.{{sample}}.cluster.log'.format(logs_d),
    conda:
        'environment.env'
    params:
        timeout = config['gurobi']['timeout'],
    threads:
        32
    shell:
        'export GRB_LICENSE_FILE={input.license}; '
        '{input.script} -s {input.segment} -o {output.cluster} -l {output.logs} -t {threads} -to {params.timeout} > {output.log}'

rule plot:
    input:
        script     = config['exec']['plot'],
        segment    = '{}/freddie.{{sample}}.segment.tsv'.format(output_d),
        cluster    = '{}/freddie.{{sample}}.cluster.tsv'.format(output_d),
        annotation = lambda wildcards: config['annotations'][config['samples'][wildcards.sample]['gtf']],
    output:
        plot_dir = directory('{}/freddie.{{sample}}.plot'.format(plots_d)),
    threads:
        32
    conda:
        'environment.env'
    shell:
        '{input.script} -a {input.annotation} -c {input.cluster} -s {input.segment} -od {output.plot_dir} -t {threads}'

rule isoforms:
    input:
        script  = config['exec']['isoforms'],
        segment = '{}/freddie.{{sample}}.segment.tsv'.format(output_d),
        cluster = '{}/freddie.{{sample}}.cluster.tsv'.format(output_d),
    output:
        isoforms = protected('{}/freddie.{{sample}}.isoforms.gtf'.format(output_d)),
    shell:
        '{input.script} -s {input.segment} -c {input.cluster} -o {output.isoforms}'

rule beds:
    input:
        script   = config['exec']['beds'],
        isoforms = '{}/freddie.{{sample}}.isoforms.gtf'.format(output_d),
    output:
        beds = directory('{}/freddie.{{sample}}.beds'.format(output_d)),
    shell:
        '{input.script} -g {input.isoforms} -od {output.beds}'

rule ref_beds:
    input:
        script_cov = config['exec']['coverage'],
        script_bed = config['exec']['beds'],
        annotation = lambda wildcards: config['annotations'][config['samples'][wildcards.sample]['gtf']],
        split      = '{}/freddie.{{sample}}.split.tsv'.format(output_d),
    output:
        ref_cov    = '{}/{{sample}}.ref_cov.tsv'.format(output_d),
        ref_beds   = directory('{}/{{sample}}.ref_beds'.format(output_d)),
    params:
        prt = lambda wildcards: config['transcript_coverage'][config['samples'][wildcards.sample]['data_type']]['per_read_threshold'],
        ptt = lambda wildcards: config['transcript_coverage'][config['samples'][wildcards.sample]['data_type']]['per_transcript_threshold'],
    shell:
        '{input.script_cov} -g {input.annotation} -s {input.split} -o {output.ref_cov} -prt {params.prt} -ptt {params.ptt}; '
        '{input.script_bed} -g {input.annotation} -c {output.ref_cov} -od {output.ref_beds}'

rule bed_intersect:
    input:
        beds     = '{}/freddie.{{sample}}.beds'.format(output_d),
        ref_beds = '{}/{{sample}}.ref_beds'.format(output_d),
    output:
        intersect = '{}/freddie.{{sample}}.intersect.tsv'.format(output_d),
    shell:
        'bedtools intersect -a <(cat {input.ref_beds}/*.bed) -b <(cat {input.beds}/*.bed) -wb'
        ' | cut -f4,8 | sort -u > "{output.intersect}"; '

rule seqpare:
    input:
        script   = config['exec']['seqpare'],
        beds     = '{}/freddie.{{sample}}.beds'.format(output_d),
        intersect = '{}/freddie.{{sample}}.intersect.tsv'.format(output_d),
        ref_beds = '{}/{{sample}}.ref_beds'.format(output_d),
    output:
        seqpare   = '{}/freddie.{{sample}}.seqpare.tsv'.format(output_d),
    shell:
        '{input.script} {input.beds} {input.ref_beds} {input.intersect} {output.seqpare}'

rule pairing:
    input:
        script   = config['exec']['pairing'],
        seqpare  = '{}/freddie.{{sample}}.seqpare.tsv'.format(output_d),
        beds     = '{}/freddie.{{sample}}.beds'.format(output_d),
        ref_beds = '{}/{{sample}}.ref_beds'.format(output_d),
    output:
        pairings = '{}/freddie.{{sample}}.pairings.tsv'.format(output_d),
    shell:
        '{input.script} -s {input.seqpare} -b {input.beds} -r {input.ref_beds} -o {output.pairings}'
