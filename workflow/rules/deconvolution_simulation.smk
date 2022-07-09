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
        simulator = lambda wildcards: str(runs.simulator[runs.index == int(wildcards.param_id)].values),
        sim_config = lambda wildcards: runs[runs.index == int(wildcards.param_id)]
    output:
        data="data/{param_id}_rep{instance_id}_data.npy",
        metadata=expand("data/{{param_id}}_rep{{instance_id}}_metadata_{model}.npy", model=config["models"])
    run:
        sim_config = params.sim_config
        alpha = np.array([float(x) for x in re.split(',',sim_config["alpha"].strip("[]").replace(" ",""))])
        sim_config["alpha"] = alpha / alpha.sum()
        sim_config["data_file"] = output.data
        sim_config["metadata_files"] = output.metadata
        sim_config["models"] = config["models"]
        name_to_simulator[sim_config["simulator"]](sim_config)

rule run_model:
    input:
        data = "data/{param_id}_rep{instance_id}_data.npy", \
        metadata = "data/{param_id}_rep{instance_id}_metadata_{model}.npy"  #lambdas
    params:
        instance="{instance_id}",
        run_config=lambda wildcards: runs[runs.index == int(wildcards.param_id)]
    output:
        "interim/{param_id}_rep{instance_id}_{model}.npy"
    run:
        model = name_to_model(params.model)
        run_config = params.run_config
        run_config["data"], run_config["metadata"] = input
        run_config["outfile"] = output[0]
        runner = model(run_config)
        runner.run_from_npy()

rule write_output:
    input:
        expand("interim/{param_id}_rep{instance_id}_{model}.npy", param_id = param_ids, instance_id = np.arange(config["reps"]), model=config["models"])
    output:
        expand("results/{name}_alpha_estimates.json", name=config["name"])
    run:
        #open each input
        for filename in input:
            alpha, i = np.load(filename, allow_pickle=True)
            with open(output[0], "a+") as outfile: #dump to file
                outfile.write(filename+"\t"+str(list(alpha))+"\t"+str(int(i))+"\n")
        #model t depth m_per_window windows_per_t thetaH thetaL lambda_high lambda_low alpha_true alpha_est
        # iterations_limit actual_iterations