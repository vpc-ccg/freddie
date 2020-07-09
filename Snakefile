configfile: 'config.yaml'

for k in config.keys():
    if not k.startswith('override_'):
        continue
    keys = k[len('override_'):].split('_')
    top_dict = eval('config{}'.format(''.join(['["{}"]'.format(x) for x in keys[:-1]])))
    assert keys[-1] in top_dict
    top_dict[keys[-1]]=config[k]

def get_abs_path(path):
    import os
    abs_path = os.popen('readlink -f {}'.format(path)).read()
    return abs_path.rstrip("\r\n")
def make_slurm():
    import os
    os.makedirs('slurm'.format(outpath), mode=0o777, exist_ok=True)

outpath = get_abs_path(config['outpath'])
make_slurm()

output_d = '{}/results'.format(outpath)
plots_d  = '{}/plots'.format(outpath)
mapped_d = '{}/mapped'.format(outpath)
logs_d   = '{}/logs'.format(outpath)

rule all:
    input:
        expand('{}/freddie.{{sample}}.sam'.format(mapped_d), sample=config['samples']),
        expand('{}/freddie.{{sample}}.split.tsv'.format(output_d), sample=config['samples']),
        expand('{}/freddie.{{sample}}.segment.tsv'.format(output_d), sample=config['samples']),
        expand('{}/freddie.{{sample}}.cluster.tsv'.format(output_d), sample=config['samples']),
        expand('{}/freddie.{{sample}}.isoforms.gtf'.format(output_d), sample=config['samples']),
        expand('{}/freddie.{{sample}}.plot/'.format(plots_d), sample=config['samples']),
        expand('{}/freddie.{{sample}}.beds/'.format(logs_d), sample=config['samples']),
        expand('{}/freddie.{{sample}}.seqpare/'.format(logs_d), sample=config['samples']),
        expand('{}/freddie.{{sample}}.similarity.txt'.format(logs_d), sample=config['samples']),

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
    conda:
        'environment.env'
    params:
        timeout = config['gurobi']['timeout'],
        logs    = '{}/freddie.{{sample}}.cluster'.format(logs_d)
    threads:
        32
    shell:
        'export GRB_LICENSE_FILE={input.license}; '
        '{input.script} -s {input.segment} -o {output.cluster} -l {params.logs} -t {threads} -to {params.timeout}'

rule isoforms:
    input:
        script  = config['exec']['isoforms'],
        cluster = '{}/freddie.{{sample}}.cluster.tsv'.format(output_d),
    output:
        isoforms = protected('{}/freddie.{{sample}}.isoforms.gtf'.format(output_d)),
    shell:
        '{input.script} -c {input.cluster} -o {output.isoforms}'

rule beds:
    input:
        script  = config['exec']['beds'],
        segment ='{}/freddie.{{sample}}.segment.tsv'.format(output_d),
        cluster = '{}/freddie.{{sample}}.cluster.tsv'.format(output_d),
        annotation = config['annotations']['gtf']
    output:
        beds_dir = directory('{}/freddie.{{sample}}.beds/'.format(logs_d)),
    shell:
        '{input.script} -a {input.annotation} -c {input.cluster} -s {input.segment} -od {output.beds_dir}'

rule seqpare:
    input:
        beds_dir = '{}/freddie.{{sample}}.beds'.format(logs_d),
    output:
        seqpare_dir = directory('{}/freddie.{{sample}}.seqpare/'.format(logs_d)),
    shell:
        """
        for t in {input.beds_dir}/*; do
            for i in $t/isos/*; do
                for j in $t/tids/*; do
                    echo -e "$(basename $i)\\t$(basename $j)";
                    paste <(seqpare $i $j|head -n1);
                done;
            done > {output.seqpare_dir}/"$(basename $t)".tsv;
        done;
        """

rule similarity:
    input:
        script  = config['exec']['similarity'],
        seqpare_dir = directory('{}/freddie.{{sample}}.seqpare/'.format(logs_d)),
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        annotation = config['annotations']['gtf'],
    output:
        similarity = protected('{}/freddie.{{sample}}.similarity.txt'.format(logs_d)),
    params:
        sim_threshold = 0.5,
        min_cov = 3,
    shell:
        '{input.script} -a {input.annotation} -r {input.reads} -sd {input.seqpare_dir} -o {output.similarity}'

rule plot:
    input:
        script = config['exec']['plot'],
        segment = '{}/freddie.{{sample}}.segment.tsv'.format(output_d),
        cluster = '{}/freddie.{{sample}}.cluster.tsv'.format(output_d),
        annotation = config['annotations']['gtf']
    output:
        plot_dir = directory('{}/freddie.{{sample}}.plot/'.format(plots_d)),
    threads:
        32
    conda:
        'environment.env'
    shell:
        '{input.script} -a {input.annotation} -c {input.cluster} -s {input.segment} -od {output.plot_dir} -t {threads}'
