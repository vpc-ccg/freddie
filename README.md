# Freddie
## Installation

The simplest way to install the dependencies is using [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/):

```bash
git clone https://baraaorabi@bitbucket.org/baraaorabi/freddie.git
cd freddie
conda env create -f environment.yml
conda activate freddie
```



## Running Freddie

There are few scripts/stages in Freddie:

- `py/freddie_align.py`: Align reads using `deSALT` mapper. (Optional since freddie accepts any SAM/BAM file for next stages)

- `py/freddie_split.py`: Partitions the reads into independent sets that can be processed in parallel
-  `py/freddie_segment.py`: Computes the canonical segmentation for each read set
-  `py/freddie_cluster.py`: Clusters the reads using their canonical segmentation representation
-  `py/freddie_isoforms.py`: Generates consensus isoforms of each cluster and outputs them in `GTF` format

### Align

```
py/freddie_align.py --reads <FASTA/FASTQ> --genome <FASTA> --out-desalt-index <DIR> --output <SAM>
```

Align takes the following arguments:

- `--reads/-r`: space separated list of FASTA/FASTQ files of the reads
- `--genome/-g`: FASTA file of the genome. Only needed to deSALT index if `--desalt-index` is not provided.

- `--out-desalt-index/-od`: Path to output directory to store deSALT index.  Only needed to deSALT index if `--desalt-index` is not provided.
- `--desalt/-d`: Path to deSALT executable. Default: `deSALT`
- `--desalt-index/-i`: Path to deSALT index if it exists.
- `--temporary-prefix/-m`: Temporary prefix for deSALT alignment files. Default: `<--output>.temp`
- `--sequencer/-s`: Sequencer option for deSALT: `null`, `ccs`, `clr`, `ont1d`, or `ont2d`. Default: `null`
- `--output/-o`: Output SAM file 

### Split (aka partition)

```
py/freddie_split.py --sam <SAM/BAM> --output <SPLIT.TSV>
```

Align takes the following arguments:

- `--sam/-s`: `SAM/BAM` file of read alignments from deSALT or any other split/splice long-read mapper.
- `--sam-format/-f`: `--sam` file format: `bam` or `sam`. Default: infers format from `--sam` file extension

- ``--output/-o`: Output TSV file of split stage. Default: `freddie_split.tsv`

### Segment

```
py/freddie_segment.py --split-tsv <SPLIT.TSV> --reads <FASTQ/FASTA> --output <SEGMENT.TSV>
```

Align takes the following arguments:

- `--split-tsv/-s`: `SPLIT.TSV` output of the split stage
- `--reads/-r`: Space separated list of FASTA/FASTQ files of the reads
- `--threads/-t`: Number of threads. Default: 1
- `--sigma/-sd`: Standard deviation parameter for the Gaussian filter. Default: 5.0
- `--threshold-rate/-tp`: Coverage percentage threshold for segments. Default: 0.90
- `--variance-factor/-vf`: Variance factor for heuristic of prefixing breakpoints. Any breakpoint with signal greater than `<-vf>` times the standard deviation plus the average of the signals will be prefixed. Default: 3.0
- `--max-problem-size/-mps`: Maximum allowed problem size after which the problem will be broken down uniformly. Default: 50
- `--min-read-support-outside`: Minimum contrasting coverage support required for a breakpoint. Default: 3
- ``--output/-o`: Output TSV file of segment stage. Default: `freddie_segment.tsv`

### Cluster

```
py/freddie_cluster.py --segment-tsv <SEGMENT.TSV> --output <CLUSTER.TSV>
```

Align takes the following arguments:

- `--segment-tsv/-s`: `SEGMENT.TSV` output of the segment stage
- `--gap-offset/-go`: Constant +/- slack used in unaligned gap condition. Default: 20
- `--epsilon/-e`: +/- ratio of length as slack used in unaligned gap condition. Default: 0.2
- `--max-rounds/-mr`: Maximum allowed number of rounds per sub-partition of a read set. Default: 30
- `--min-isoform-size/-is`: Minimum read support allowed for an isoform. Default: 3
- `--timeout/-to`: Gurobi timeout per isoform in minutes. Default: 4
- `--threads/-t`: Number of threads. Default: 1
- `--logs-dir/-l`: Directory path where logs will be outputted. Default: No logging
- `--output/-o`: Output TSV file of cluster stage. Default: `freddie_cluster.tsv`

### Isoforms

```
py/freddie_isoforms.py --segment-tsv <SEGMENT.TSV> --cluster-tsv <CLUSTER.TSV> --output <ISOFORMS.GTF>
```

Align takes the following arguments:

- `--segment-tsv/-s`: `SEGMENT.TSV` output of the segment stage
- `--cluster-tsv/-s`: `CLUSTER.TSV` output of the cluster stage
- `--output/-o`: Output GTF file of isoforms stage. Default: `freddie_isoforms.gtf`


