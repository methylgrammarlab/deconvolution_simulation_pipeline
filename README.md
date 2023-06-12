# Simulation pipeline for deconvolution models


The pipeline uses two types of configuration files:
1) config.json - general parameters such as run name, models to use, number of replicates
2) tabular config - each row represents an in-silico mixture. 
The tabular config files used to create paper figures can be found under "resources".

In-house packages required to run this pipeline: 
1) [epiread-tools](https://github.com/methylgrammarlab/epiread-tools) 
2) [deconvolution_models](https://github.com/methylgrammarlab/deconvolution_models) (should automatically install epiread-tools)

## Usage

Clone the repository, adjust config.json as necessary. From the project root simply run
```
snakemake --cores 1
```
Output files will be created under "results". 