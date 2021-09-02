# Freddie Benchmarking
The branch is to be used for regenerating the benchmarking results of Freddie's paper.

## Installation
To start off, please install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Then, download the repo and switch to `benchmarking` branch:
```bash
(base)$ git clone https://github.com/vpc-ccg/freddie.git freddie
(base)$ cd freddie
(base)$ git checkout benchmarking
(base)$ git submodule update --init --recursive
```

## Gurobi license
The clustering stage of Freddie uses Gurobi solver which needs a license to use.
If your affliation is academic, you can cost-free obtain a license [here](https://www.gurobi.com/downloads/end-user-license-agreement-academic/).
Make sure to update the license path in `config.yaml` to point to the installed license file.


## Generating the simulated reads:
When downloading the repo, you will have downloaded the `LTR-sim` repo too. 
We will use this to generate the simulated reads.

First off, we will go the the `LTR-sim` submodule and install its Conda environment:
```bash
(base)$ cd extern/LTR-sim/
(base)$ conda env create -f env.yml
(base)$ conda activate ltr-sim
```

Here, you can edit the `config.yaml` (in `LTR-sim` folder) including editing the samples and their simulated properties.
Note that each sample depends on using a pre-existing real long-reads dataset.
A real dataset can either be provided directly or it can be downloaded using its SRA number.
Also, you may change the `batches` to speedup the parallalization of Badread simulator.

To make things easier, the transcript read expression counts for two datasets have already been added to the repo so you don't have to download them:
- 22Rv1 cell line (SRA #SRR14374285): `output/train/22Rv1.expression_rate.tsv`
- Fruit fly (SRA #ERR3588905): `output/train/ERR3588905.expression_rate.tsv`

You will not need to download the real raw reads for these two datasets to reproduce the simulated results in the paper.

Once you configured simulation, run Snakemake to generate the simulation reads:
```bash
(ltr-sim)$ snakemake -j <threads>
(ltr-sim)$ conda deactivate
(base)$ cd ../..
```

Note that you can change [Snakemake's execution parameters](https://snakemake.readthedocs.io/en/stable/executing/cli.html) to fit your server properties.

### Running the isoform detection tools
The config file already points to the two simulated datasets.
If you added any other datasets, or if you have a real dataset that you want to test, you need to add them to the `config.yaml` file that is in the root directory of the repo.


Now, first thing you need to do is to create an environment for running the tools:
```bash
(base)$ conda env create -f envs/run_tools.yml
(base)$ conda activate freddie_bench_run_tools
```

Then create the environments for each specific tool:
```bash
(freddie_bench_run_tools)$ snakemake -s Snakefile-run_tools --use-conda --conda-create-envs-only -j  <threads>
```


Now, you can run the tools:
```bash
(freddie_bench_run_tools)$ snakemake -s Snakefile-run_tools --use-conda --conda-create-envs-only -j  <threads>
```

And you can run the accuracy assessment scripts:
```bash
(freddie_bench_run_tools)$ snakemake -s Snakefile-accuracy -j <threads>
```
