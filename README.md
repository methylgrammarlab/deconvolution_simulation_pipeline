# Simulation pipeline for deconvolution models

The pipeline uses two types of configuration files:
1) config.json - general parameters such as run name, models to use, number of replicates
2) tabular config - each row represents a simulation. It includes which simulator to use and all
simulator input parameters. The tabular config files used to create paper figures can be found under "resources".

Two in-house packages are required to run this pipeline: 
1) https://github.com/methylgrammarlab/epiread-tools (install this first)
2) https://github.com/methylgrammarlab/deconvolution_models 

## Usage

Clone the repository, adjust config.json as necessary. From the project root simply run
```
snakemake --cores 1
```
Output files will be created under "results". 