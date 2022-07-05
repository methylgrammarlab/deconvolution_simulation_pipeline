from deconvolution_models.main import Celfie, CelfiePlus, Epistate, EpistatePlus
import numpy as np
import pandas as pd

runs = pd.read_csv(config["run_file"]).set_index("param_id")
param_ids = runs.index
name_to_model = {"celfie":Celfie, "celfie-plus":CelfiePlus, "epistate":Epistate, "epistate-plus":EpistatePlus}


rule simulate_data: #should I separate atlas from mixture?

    params:
        instance = "{instance_id}",
        simulator = lambda wildcards: str(runs.simulator[runs.index == int(wildcards.param_id)].values),
        sim_config = lambda wildcards: runs[runs.index == int(wildcards.param_id)]
    output:
        data="data/{param_id}_rep{instance_id}_data.npy",
        metadata="data/{param_id}_rep{instance_id}_metadata.npy"
    shell:
        """python3 epistate_simulator.py {params.t} {params.depth} {params.m_per_window} {params.windows_per_t} """ +\
            """{params.thetaL} {params.thetaH} {params.lambda_high} {params.lambda_low} """+\
                """{output.data} {output.metadata} {params.alpha} 2> {log} """


rule run_model:
    input:
        data = "data/{param_id}_rep{instance_id}_data.npy", \
        metadata = "data/{param_id}_rep{instance_id}_metadata.npy"  #lambdas
    params:
        instance="{instance_id}",
        model="{model}",
        run_config=lambda wildcards: runs[runs.index == int(wildcards.param_id)]
    output:
        "interim/{param_id}_rep{instance_id}_{model}.npy"
    run:
        model = name_to_model(params.model)
        runner = model(params.run_config) #TODO add input and output files to config
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