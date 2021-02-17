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

def get_which(command):
    import subprocess
    result = subprocess.run(['which', command], stdout=subprocess.PIPE)
    return result.stdout.decode('ascii').rstrip()

outpath      = config['outpath'].rstrip('/')
preprocess_d = '{}/preprocess'.format(outpath)
workspace_d  = '{}/workspace'.format(outpath)
output_d     = '{}/output'.format(outpath)

gnu_time = get_which('time')
bam2Bed12 = '{}/bin/bam2Bed12.py'.format(get_which('flair.py').rstrip('/flair.py'))

make_slurm()

mappers = [
    'desalt',
    'minimap2'
]
tools = [
    'freddie',
    'stringtie',
    'flair',
]

gtf_sample_rates = [
    1.00,
    0.75,
    0.50,
    0.25,
    0.01,
]
tools_stats = list()
for mapper in mappers:
    for tool in tools:
        if tool == 'flair':
            for r in gtf_sample_rates:
                tools_stats.append('flair.{:.2f}.{}'.format(r,mapper))
        else:
            tools_stats.append('{}.{}'.format(tool, mapper))

stats_outfiles = [
    'isoforms.gtf',
]

rule all:
    input:
        [config['references'][s]['desalt_idx'] for s in config['references']],
        expand('{}/{{sample}}.{{mapper}}.{{extension}}'.format(preprocess_d),
            sample=config['samples'],
            mapper=mappers,
            extension=['sam','bam','bam.bai']
        ),
        
        # expand('{}/{{sample}}/{{tool}}.{{mapper}}.{{extension}}'.format(output_d),
        #     tool=tools,
        #     mapper=mappers,
        #     sample=config['samples'],
        #     extension=['time.tsv',]),
        # expand('{}/{{sample}}/{{tool}}.{{extension}}'.format(output_d),
        #     tool=tools_stats,
        #     sample=config['samples'],
        #     extension=stats_outfiles),
        # expand('{}/{{sample}}/seqpare_stats.tsv'.format(output_d), sample=config['samples']),

rule minimap2:
    input:
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        genome = lambda wildcards: config['references'][config['samples'][wildcards.sample]['ref']]['genome'],
    output:
        sam=protected('{}/{{sample}}.minimap2.sam'.format(preprocess_d)),
        bam=protected('{}/{{sample}}.minimap2.bam'.format(preprocess_d)),
        bai=protected('{}/{{sample}}.minimap2.bam.bai'.format(preprocess_d)),
    threads:
        32
    shell:
        'minimap2 -a -x splice -t {threads} {input.genome} {input.reads} > {output.sam} && '
        '  samtools sort -T {output.bam}.tmp -m 2G -@ {threads} -O bam {output.sam} > {output.bam} && '
        '  samtools index {output.bam} '

rule desalt_index:
    input:
        genome = lambda wildcards: config['references'][wildcards.species]['genome'],
    output:
        index  = directory('test/mapping/{species}.dna.desalt_idx'),
    wildcard_constraints:
        species = '|'.join(s for s in config['references']),
    shell:
        'deSALT index {input.genome} {output.index}'

rule desalt:
    input:
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        index  = lambda wildcards: config['references'][config['samples'][wildcards.sample]['ref']]['desalt_idx'],
    output:
        sam=protected('{}/{{sample}}.desalt.sam'.format(preprocess_d)),
        bam=protected('{}/{{sample}}.desalt.bam'.format(preprocess_d)),
        bai=protected('{}/{{sample}}.desalt.bam.bai'.format(preprocess_d)),
    params:
        seq_type = lambda wildcards: config['samples'][wildcards.sample]['seq_type'],
    threads:
        32
    shell:
        'py/freddie_align.py -r {input.reads} -i {input.index}/ -s {params.seq_type} -o {output.sam} -t {threads} &&'
        '  samtools sort -T {output.bam}.tmp -m 2G -@ {threads} -O bam {output.sam} > {output.bam} && '
        '  samtools index {output.bam} '

rule freddie:
    input:
        sam   = '{}/{{sample}}.{{mapper}}.sam'.format(preprocess_d),
        reads =  lambda wildcards: config['samples'][wildcards.sample]['reads'],
    output:
        split   = protected('{}/{{sample}}/freddie.{{mapper}}/split.tsv'.format(workspace_d)),
        segment = protected('{}/{{sample}}/freddie.{{mapper}}/segment.tsv'.format(workspace_d)),
        cluster = protected('{}/{{sample}}/freddie.{{mapper}}/cluster.tsv'.format(workspace_d)),
    params:
        gnu_time     = gnu_time,
        gnu_time_fmt = '"%e\\t%U\\t%M"',
        gurobi       = config['gurobi']['license']
    wildcard_constraints:
        mapper='desalt|minimap2'
    shell:
        '{params.gnu_time} -f {params.gnu_time_fmt} -o {output.time} bash -c "'
        ' py/freddie_split.py -s {input.sam} -o {output.split}'
        '" && '
        '{params.gnu_time} -f {params.gnu_time_fmt} -a -o {output.time} bash -c "'
        ' py/freddie_segment.py -s {output.split} -r {input.reads} -o {output.segment} '
        '" && '
        '{params.gnu_time} -f {params.gnu_time_fmt} -a -o {output.time} bash -c "'
        ' export GRB_LICENSE_FILE={params.gurobi} && '
        ' py/freddie_cluster.py -s {output.segment} -o {output.cluster} -to 15 > {output.cluster}.log'
        '"'

rule freddie_gtf:
    input:
        segment = protected('{}/{{sample}}/freddie.{{mapper}}/segment.tsv'.format(workspace_d)),
        cluster = protected('{}/{{sample}}/freddie.{{mapper}}/cluster.tsv'.format(workspace_d)),
    output:
        gtf     = protected('{}/{{sample}}/freddie.{{mapper}}.isoforms.gtf'.format(output_d)),
    wildcard_constraints:
        mapper='desalt|minimap2'
    shell:
        'py/freddie_isoforms.py -s {input.segment} -c {input.cluster} -o {output.gtf}'

rule stringtie:
    input:
        bam = '{}/{{sample}}.{{mapper}}.bam'.format(preprocess_d),
    output:
        gtf  = protected('{}/{{sample}}/stringtie.{{mapper}}.isoforms.gtf'.format(output_d)),
        time = protected('{}/{{sample}}/stringtie.{{mapper}}.time.tsv'.format(output_d)),
    wildcard_constraints:
        mapper='desalt|minimap2'
    params:
        gnu_time     = gnu_time,
        gnu_time_fmt = '"%e\\t%U\\t%M"',
    shell:
        '{params.gnu_time} -f {params.gnu_time_fmt} -o {output.time} '
        ' stringtie -p 1 -L -o {output.gtf} {input.bam}'

rule sampled_gtf:
    input:
        gtf = lambda wildcards: config['annotations'][config['samples'][wildcards.sample]['gtf']],
    output:
        gtf ='{}/{{sample}}.sampled_at.{{gtf_sample_rate}}.gtf'.format(preprocess_d),
    run:
        import numpy as np
        print(input.gtf)
        np.random.seed(42)
        isoform_ids = set()
        isoform_id_gid = dict()
        for l in open(input.gtf):
            if l[0]=='#':
                continue
            l = l.rstrip().split('\t')
            if l[2]!='transcript':
                continue
            info = l[8]
            info = [x.strip().split(' ') for x in info.strip(';').split(';')]
            info = {x[0]:x[1].strip('"') for x in info}
            isoform_ids.add(info['transcript_id'])
            isoform_id_gid[info['transcript_id']] = info['gene_id']
        print(len(isoform_ids))
        isoform_ids = set(np.random.choice(
            list(isoform_ids),
            int(len(isoform_ids)*float(wildcards.gtf_sample_rate)),
            replace=False,
        ))
        print(len(isoform_ids))
        gene_ids = {isoform_id_gid[isoform_id] for isoform_id in isoform_ids}
        print(len(gene_ids))
        out_file = open(output.gtf, 'w+')
        for line in open(input.gtf):
            if line[0]=='#':
                out_file.write(line)
                continue
            l = line.rstrip().split('\t')
            info = l[8]
            info = [x.strip().split(' ') for x in info.strip(';').split(';')]
            info = {x[0]:x[1].strip('"') for x in info}
            if l[2]=='gene' and info['gene_id'] in gene_ids:
                out_file.write(line)
            elif l[2]!='gene' and info['transcript_id'] in isoform_ids:
                out_file.write(line)
        out_file.close()

rule flair:
    input:
        bam    ='{}/{{sample}}.{{mapper}}.bam'.format(preprocess_d),
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        genome = lambda wildcards: config['references'][config['samples'][wildcards.sample]['gtf']],
        gtf    ='{}/{{sample}}.sampled_at.{{gtf_sample_rate}}.gtf'.format(preprocess_d),
    output:
        p_bed  = protected('{}/{{sample}}/flair.{{gtf_sample_rate}}.{{mapper}}/input.bed'.format(workspace_d)),
        c_bed  = protected('{}/{{sample}}/flair.{{gtf_sample_rate}}.{{mapper}}/correct_all_corrected.bed'.format(workspace_d)),
        i_bed  = protected('{}/{{sample}}/flair.{{gtf_sample_rate}}.{{mapper}}/correct_all_inconsistent.bed'.format(workspace_d)),
        t_bed  = protected('{}/{{sample}}/flair.{{gtf_sample_rate}}.{{mapper}}/collapse.isoforms.bed'.format(workspace_d)),
        fasta  = protected('{}/{{sample}}/flair.{{gtf_sample_rate}}.{{mapper}}/collapse.isoforms.fa'.format(workspace_d)),
        gtf  = protected('{}/{{sample}}/flair.{{gtf_sample_rate}}.{{mapper}}.isoforms.gtf'.format(output_d)),
        time = protected('{}/{{sample}}/flair.{{gtf_sample_rate}}.{{mapper}}.time.tsv'.format(output_d)),
    params:
        gnu_time        = gnu_time,
        gnu_time_fmt    = '"%e\\t%U\\t%M"',
        bam2Bed12       = bam2Bed12,
        correct_prefix  = '{}/{{sample}}/flair.{{gtf_sample_rate}}.{{mapper}}/correct'.format(workspace_d),
        collapse_prefix = '{}/{{sample}}/flair.{{gtf_sample_rate}}.{{mapper}}/collapse'.format(workspace_d),
    wildcard_constraints:
        mapper='desalt|minimap2'
    shell:
        ' {params.bam2Bed12} -i {input.bam} > {output.p_bed}  && '
        '{params.gnu_time} -f {params.gnu_time_fmt} -o {output.time} bash -c "'
        ' flair.py correct  -t 1 -q {output.p_bed} -o {params.correct_prefix} -g {input.genome} -f {input.gtf} &&'
        ' flair.py collapse -t 1 -q {output.c_bed} -o {params.collapse_prefix} -g {input.genome} -f {input.gtf} -r {input.reads} '
        '"; '
        'mv {params.collapse_prefix}.isoforms.gtf {output.gtf}'

rule ref_beds:
    input:
        sam  = '{}/{{sample}}.{{mapper}}.sam'.format(preprocess_d),
        script_cov = config['exec']['coverage'],
        script_bed = config['exec']['beds'],
        annotation = lambda wildcards: config['annotations'][config['samples'][wildcards.sample]['gtf']],
    output:
        split      = '{}/{{sample}}.{{mapper}}.split.tsv'.format(preprocess_d),
        ref_cov    = '{}/{{sample}}.{{mapper}}.ref_cov.tsv'.format(preprocess_d),
        ref_beds   = directory('{}/{{sample}}.{{mapper}}.ref_beds'.format(preprocess_d)),
    params:
        prt = 0.0,
        ptt = 3,
    wildcard_constraints:
        mapper = 'desalt|minimap2'
    shell:
        'py/freddie_split.py -s {input.sam} -o {output.split} && '
        '{input.script_cov} -g {input.annotation} -s {output.split} -o {output.ref_cov} -prt {params.prt} -ptt {params.ptt} && '
        '{input.script_bed} -g {input.annotation} -c {output.ref_cov} -od {output.ref_beds}'

rule tool_beds:
    input:
        script   = config['exec']['beds'],
        isoforms = '{}/{{sample}}/{{tool}}.{{mapper}}.isoforms.gtf'.format(output_d),
    output:
        tool_beds = directory('{}/{{sample}}/{{tool}}.{{mapper}}.isoforms.beds'.format(workspace_d)),
    shell:
        '{input.script} -g {input.isoforms} -od {output.tool_beds}'


rule bed_intersect:
    input:
        ref_beds  = '{}/{{sample}}.{{mapper}}.ref_beds'.format(preprocess_d),
        tool_beds = '{}/{{sample}}/{{tool}}.{{mapper}}.isoforms.beds'.format(workspace_d),
    output:
        intersect = '{}/{{sample}}/{{tool}}.{{mapper}}.intersect.tsv'.format(output_d),
    shell:
        'bedtools intersect -a <(cat {input.ref_beds}/*.bed) -b <(cat {input.tool_beds}/*.bed) -wb'
        ' | cut -f4,8 | sort -u > "{output.intersect}"; '

rule seqpare:
    input:
        script    = config['exec']['seqpare'],
        intersect = '{}/{{sample}}/{{tool}}.{{mapper}}.intersect.tsv'.format(output_d),
        ref_beds  = '{}/{{sample}}.{{mapper}}.ref_beds'.format(preprocess_d),
        tool_beds = '{}/{{sample}}/{{tool}}.{{mapper}}.isoforms.beds'.format(workspace_d),
    output:
        seqpare   = '{}/{{sample}}/{{tool}}.{{mapper}}.seqpare.tsv'.format(output_d),
    shell:
        '{input.script} {input.tool_beds} {input.ref_beds} {input.intersect} {output.seqpare}'

rule pairing:
    input:
        script    = config['exec']['pairing'],
        seqpare   = '{}/{{sample}}/{{tool}}.{{mapper}}.seqpare.tsv'.format(output_d),
        tool_beds = '{}/{{sample}}/{{tool}}.{{mapper}}.isoforms.beds'.format(workspace_d),
        ref_beds  = '{}/{{sample}}.{{mapper}}.ref_beds'.format(preprocess_d),
    output:
        pairings = '{}/{{sample}}/{{tool}}.{{mapper}}.pairings.tsv'.format(output_d),
    shell:
        '{input.script} -s {input.seqpare} -b {input.tool_beds} -r {input.ref_beds} -o {output.pairings} -t 1'

rule stats:
    input:
        script   = config['exec']['stats'],
        pairing_files = ['{}/{{sample}}/{}.pairings.tsv'.format(output_d, tool) for tool in tools_stats],
    output:
        stats = '{}/{{sample}}/seqpare_stats.tsv'.format(output_d),
    shell:
        '{input.script} {input.pairing_files} > {output.stats}'
