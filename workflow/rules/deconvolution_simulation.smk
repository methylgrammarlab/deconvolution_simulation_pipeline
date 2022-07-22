import sys
sys.path.append("/Users/ireneu/PycharmProjects/deconvolution_models")
from deconvolution_models.main import Celfie, CelfiePlus, Epistate, EpistatePlus
from scripts.epistate_simulator import main as epistate_simulator
from scripts.random_simulator import main as random_simulator
import numpy as np
import pandas as pd
import re

runs = pd.read_csv(config["run_file"]).set_index("param_id")
param_ids = runs.index
name_to_model = {"celfie":Celfie, "celfie-plus":CelfiePlus, "epistate":Epistate, "epistate-plus":EpistatePlus}
name_to_simulator = {"epistate":epistate_simulator, "random":random_simulator}

rule simulate_data: #should I separate atlas from mixture?
    params:
        instance = "{instance_id}",
        simulator = lambda wildcards: str(runs.simulator[runs.index == int(wildcards.param_id)].values[0]),
        sim_config = lambda wildcards: runs[runs.index == int(wildcards.param_id)].iloc[0,:].to_dict()
    output:
        data=expand("data/{name}/{{param_id}}_rep{{instance_id}}_data.npy", name=config["name"]),
        metadata=expand("data/{name}/{{param_id}}_rep{{instance_id}}_metadata_{model}.npy",
            model=config["models"], name=config["name"])

    run:
        sim_config = params.sim_config
        alpha = np.array([float(x) for x in re.split(',',sim_config["true_alpha"].strip("[]").replace(" ",""))])
        sim_config["true_alpha"] = alpha / alpha.sum()
        sim_config["data_file"] = output.data[0]
        sim_config["metadata_files"] = output.metadata
        sim_config["models"] = config["models"]
        name_to_simulator[params.simulator](sim_config)

rule run_model:
    input:
        data=expand("data/{name}/{{param_id}}_rep{{instance_id}}_data.npy", name=config["name"]),
        metadata=expand("data/{name}/{{param_id}}_rep{{instance_id}}_metadata_{{model}}.npy" ,name=config["name"])
    params:
        instance="{instance_id}",
        run_config=lambda wildcards: runs[runs.index == int(wildcards.param_id)].iloc[0,:].to_dict()
    output:
        temp("interim/{param_id}_rep{instance_id}_{model}.npy")
    run:
        model = name_to_model[wildcards.model]
        run_config = params.run_config
        run_config["data_file"], run_config["metadata_file"] = input.data[0], input.metadata[0]
        run_config["outfile"] = output[0]
        runner = model(run_config)
        runner.run_from_npy()

rule write_output:
    input:
        expand("interim/{param_id}_rep{instance_id}_{model}.npy", param_id = param_ids,
            instance_id = np.arange(config["reps"]), model=config["models"])
    output:
        expand("results/{name}_alpha_estimates.tsv", name=config["name"])
    run:
        #open each input
        for filename in input:
            alpha, i = np.load(filename, allow_pickle=True)
            with open(output[0], "a+") as outfile: #dump to file
                outfile.write(filename+"\t"+str(list(alpha))+"\t"+str(int(i))+"\n")
        #model t depth m_per_window windows_per_t thetaH thetaL lambda_high lambda_low alpha_true alpha_est
        # iterations_limit actual_iterations