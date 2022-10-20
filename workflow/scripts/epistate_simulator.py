# ###################################################
#
# MIT License
#
# Copyright (c) 2022 irene unterman and ben berman
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ###################################################

import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
from epiread_tools.naming_conventions import *
import numpy as np
import re
from itertools import cycle


'''
This simulator generates reads according to the two-epistate model.
We have a list of regions. Each region contains m CpGs.
Each region has two epistates associated with it - theta_high and theta_low.
Each epistate is a series of methylation probabilities per CpG.
Each methylation read can only originate from one epistate. 
The probability of a read in a given cell type t, belonging to theta_high is lambda_t.
'''

#%%

def adjust_length(var, req_length):
    if isinstance(var, float) or isinstance(var, int): #repeat if necessary.
        return np.array([var]*req_length)
    elif len(var)==req_length:
        return np.array(var)
    else:
        raise ValueError("var not num nor in required length")


def generate_atlas(config): #get lambdas, thetas
    '''
    an atlas is defined by the epistates for each region
    anf the probability of the epistates per cell type
    :param config: params for atlas generation
    :return: theta, lambdas
    '''
    theta_high = []
    theta_low = []
    lambda_t = []
    high = adjust_length(config["theta_high"], config["m_per_region"])
    low = adjust_length(config["theta_low"], config["m_per_region"])

    l_t = np.ones((config["t"]))
    l_t.fill(config["lambda_high"])
    for t in range(config["t"]):
        for i in range(config["regions_per_t"]):
            theta_high.append(high)
            theta_low.append(low)
            l = l_t.copy()
            l[t] = config["lambda_low"]
            lambda_t.append(l)
    return theta_high, theta_low, lambda_t

#%%

def generate_mixture(atlas, alpha, coverage): #TODO: re-implement random regions
    '''

    :param atlas: theta_high, theta_low, lambda_t
    :param alpha: proportions in mixture
    :param coverage: coverage of mixture
    :return: mixture
    '''
    theta_high, theta_low, lambda_t = atlas
    t, m_per_region = lambda_t[0].shape[0], theta_low[0].shape[0]
    n_regions = len(lambda_t)
    mixture = []

    for i in range(n_regions):
        mix = np.zeros((coverage,m_per_region))
        true_z = np.random.choice(np.arange(t), size=coverage, p=alpha, replace=True)
        for j in range(coverage):
            mu_i = np.random.binomial(1, p=lambda_t[i][true_z[j]])
            theta = (theta_low, theta_high)[mu_i][i]
            row = np.array([UNMETHYLATED, METHYLATED])[np.random.binomial(n=1, p=theta, size=m_per_region)]
            mix[j, :] = row
        mixture.append(mix)
    mixture = [np.vstack(x) if len(x) else np.array([]) for x in mixture]
    return mixture

def atlas_to_beta_matrices(atlas):
    thetaH, thetaL, lambdas = atlas
    t, m_per_region = lambdas[0].shape[0], thetaL[0].shape[0]
    beta_tm = []
    for i in range(len(lambdas)):
        beta = lambdas[i][:, np.newaxis]*thetaH[i]+ (1-lambdas[i])[:, np.newaxis]*thetaL[i]
        beta_tm.append(beta)
    return beta_tm

def save_mixture(data_file, reads):
    np.save(data_file, reads, allow_pickle=True)

def write_celfie_output(output_file, beta_tm, atlas_coverage=1000): #no summing
    y = np.hstack([(x*atlas_coverage) for x in beta_tm])
    y_depths = np.ones((y.shape))
    y_depths.fill(atlas_coverage)
    np.save(output_file, [y, y_depths], allow_pickle=True)

def write_celfie_plus_output(output_file, beta_tm):
    np.save(output_file, beta_tm, allow_pickle=True)

def write_epistate_output(output_file, atlas):
    thetaH, thetaL, lambdas = atlas
    np.save(output_file, [thetaH, thetaL, lambdas], allow_pickle=True)

def main(config):
    atlas = generate_atlas(config)
    beta = atlas_to_beta_matrices(atlas)
    reads = generate_mixture(atlas, config["true_alpha"], config["coverage"])
    save_mixture(config["data_file"], reads)
    for model, outfile in zip(config["models"], config["metadata_files"]):
        if model == "celfie":
            write_celfie_output(outfile, beta)
        elif model == "celfie-plus":
            write_celfie_plus_output(outfile, beta)
        else:
            write_epistate_output(outfile, atlas)

# config = {'simulator': 'epistate', 't': 25, 'coverage': 10, 'm_per_region': 6, 'regions_per_t': 40, 'num_iterations': 1000,
#           'atlas_coverage': 1000, 'random_restarts': 1, "stop_criterion":0.001,
#           'theta_high': 0.66, 'theta_low': 0.34, 'lambda_high': 1, 'lambda_low': 0,
#           "models":["celfie-plus", "celfie", "epistate-plus","epistate"],
#           "data_file": "/Users/ireneu/PycharmProjects/deconvolution_simulation_pipeline/data/length_variation/16_rep0_data.npy",
#           "metadata_files": ["/Users/ireneu/PycharmProjects/deconvolution_simulation_pipeline/data/16_rep0_metadata_celfie-plus.npy",
# "/Users/ireneu/PycharmProjects/deconvolution_simulation_pipeline/data/16_rep0_metadata_celfie.npy",
# "/Users/ireneu/PycharmProjects/deconvolution_simulation_pipeline/data/16_rep0_metadata_epistate-plus.npy",
# "/Users/ireneu/PycharmProjects/deconvolution_simulation_pipeline/data/16_rep0_metadata_epistate.npy"],
#
#           "outfile":None}
#
#
# a = np.arange(1,26)
# b = a/np.sum(a)
# config["true_alpha"] = b
# main(config)