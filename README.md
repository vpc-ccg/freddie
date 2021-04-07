# Freddie

## Running Freddie using Snakemake

The whole Freddie pipeline is readily available using Snakemake.
You can check the pipeile at the `Snakefile` and you can add samples and other paths (e.g. Gurobi licence, reference genome) to run Freddie on in the `config.yaml` file.
After editing `config.yaml`, you can run Snakemake with your specific settings, just make sure to use `--use-conda` to have all the requirements installed on the fly.
Note that the cluster stage uses Gurobi solver which needs a license to use.
If your affliation is academic, you can cost-free obtain a license [here](https://www.gurobi.com/downloads/end-user-license-agreement-academic/).
Make sure to update the license path in `config.yaml` to point to the installed license file.


## Running Freddie manually

### Installation

The simplest way to install the dependencies is using [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/):

```
git clone https://github.com/vpc-ccg/freddie.git
cd freddie
conda env create -f envs/freddie.yml
conda activate freddie
```

There are few scripts/stages in Freddie:

- `py/freddie_split.py`: Partitions the reads into independent sets that can be processed in parallel
- `py/freddie_segment.py`: Computes the canonical segmentation for each read set
- `py/freddie_cluster.py`: Clusters the reads using their canonical segmentation representation
- `py/freddie_isoforms.py`: Generates consensus isoforms of each cluster and outputs them in `GTF` format

### Align

```
minimap2 -a -x splice -t {threads} {genome FASTA} {reads FASTA/FASTQ} > {SAM}
```


### Sort
Before running split stage, the SAM file needs to be sorted and indexed using SAMtools

```
samtools sort {SAM} -m {memory per thread e.g. 2GB} -@ {threads} -O bam > {BAM}
samtools index {BAM}
```

### Split (aka partition)

```
py/freddie_split.py --bam {sorted BAM} --outdir {SPLIT} -t {threads}
```

Align takes the following arguments:

- `--reads/-r`: Space separated paths to reads in FASTQ or FASTA
- `--bam/-b`: `BAM` file of read alignments from a split/splice long-read mapper that are position sorted and indexed.
- `--outdir/-o`: Output TSV file of split stage. Default: `freddie_split/`

### Segment

```
py/freddie_segment.py -s {SPLIT} --outdir  {SEGMENT} -t {threads}
```

Align takes the following arguments:

- `--split-dir/-s`: `SPLIT` output directory of the split stage
- `--threads/-t`: Number of threads. Default: 1
- `--sigma/-sd`: Standard deviation parameter for the Gaussian filter. Default: 5.0
- `--threshold-rate/-tp`: Coverage percentage threshold for segments. Default: 0.90
- `--variance-factor/-vf`: Variance factor for heuristic of prefixing breakpoints. Any breakpoint with signal greater than `-vf` times the standard deviation plus the average of the signals will be prefixed. Default: 3.0
- `--max-problem-size/-mps`: Maximum allowed problem size after which the problem will be broken down uniformly. Default: 50
- `--min-read-support-outside`: Minimum contrasting coverage support required for a breakpoint. Default: 3
- ``--outdir/-o`: Output directory of segment stage. Default: `freddie_segment/`

### Cluster
The cluster stage uses Gurobi solver which needs a license to use.
If your affliation is academic, you can cost-free obtain a license [here](https://www.gurobi.com/downloads/end-user-license-agreement-academic/).


```
export GRB_LICENSE_FILE={path to Gurobi v9 license}
py/freddie_cluster.py --segment-dir {SEGMENT} --outdir {CLUSTER}
```

Align takes the following arguments:

- `--segment-dir/-s`: `SEGMENT` output directory of the segment stage
- `--gap-offset/-go`: Constant +/- slack used in unaligned gap condition. Default: 20
- `--epsilon/-e`: +/- ratio of length as slack used in unaligned gap condition. Default: 0.2
- `--max-rounds/-mr`: Maximum allowed number of rounds per sub-partition of a read set. Default: 30
- `--min-isoform-size/-is`: Minimum read support allowed for an isoform. Default: 3
- `--timeout/-to`: Gurobi timeout per isoform in minutes. Default: 4
- `--threads/-t`: Number of threads. Default: 1
- `--logs-dir/-l`: Directory path where logs will be outputted. Default: No logging
- `--outdir/-o`: Output directory of cluster stage. Default: `freddie_cluster/`

### Isoforms

```
py/freddie_isoforms.py --segment-dir {SEGMENT} --cluster-dir {CLUSTER} --output {ISOFORMS.GTF} -t {threads}
```

Align takes the following arguments:

- `--segment-dir/-s`: `SEGMENT` output directory of the segment stage
- `--cluster-dir/-s`: `CLUSTER` output directory of the cluster stage
- `--output/-o`: Output GTF file of isoforms stage. Default: `freddie_isoforms.gtf`
- `--threads/-t`: Number of threads. Default: 1


