#################################################
# FILE: random_simulator.py
# WRITER: Irene Unterman
# DESCRIPTION: simulate methylation sequencing
# with an atlas of random methylation probabilities
#
#################################################
import numpy as np
import sys
sys.path.append("/Users/ireneu/PycharmProjects/epiread-tools")
from epiread_tools.naming_conventions import *


def generate_atlas(t, m_per_window, n_windows):
    '''
    simulate pure cell types
    :param t: number of cell types to generate
    :param m_per_window: number of CpGs in each region
    :param n_windows: number of regions
    :return: list of beta value arrays
    '''
    atlas = []
    for i in range(n_windows):
        atlas.append(np.random.uniform(
            0, 1, size=(t, m_per_window)
        ))
    return atlas

def generate_mixture(atlas, alpha, coverage):
    '''
    sample WGBS sequences from an atlas of methylation probabilities
    :param atlas: t cell types by m CpGs. filled with beta values
    :param alpha: proportions of cell types to sample. ordered by atlas
    :param coverage: amount of reads to sample, each will be the entire
                    length of the atlas (m CpGs)
    :return: depthXm methylayion np array
    '''
    t, m = atlas[0].shape #cell types, cpgs

    n_regions = len(atlas)
    n_reads = coverage*n_regions

    mixture = [[] for x in range(n_regions)]
    region_indicator = np.random.choice(np.arange(n_regions), size=n_reads, replace=True)
    true_z = np.random.choice(np.arange(t), size=n_reads, p=alpha, replace=True)

    for i in range(n_reads):
        beta = atlas[region_indicator[i]][true_z[i],:]
        row = np.array([UNMETHYLATED, METHYLATED])[np.random.binomial(n=1, p=beta, size=m)]
        mixture[region_indicator[i]].append(row)
    mixture = [np.vstack(x) if len(x) else np.array([]) for x in mixture]

    return mixture

def make_thetas_and_lambdas(atlas, t, windows_per_t, Lhigh=1, Llow=0):
    thetaH = []
    thetaL = []
    lambdas = []
    for cell_type in range(t):
        lambda_t = np.ones((t))
        lambda_t.fill(Lhigh)
        lambda_t[cell_type] = Llow
        lambdas.extend([lambda_t]*windows_per_t)
    for i in range(len(atlas)):
        thetaH.append((atlas[i][lambdas[i]==Lhigh]).mean(axis=0))
        thetaL.append((atlas[i][lambdas[i]==Llow]).mean(axis=0))
    return thetaH, thetaL, lambdas


def generate_data(config):
    atlas = generate_atlas(config["t"], config["m_per_region"], config["t"]*config["regions_per_t"])
    reads = generate_mixture(atlas, config["true_alpha"], config["coverage"])
    thetaH, thetaL, lambdas =  make_thetas_and_lambdas(atlas, config["t"], config["regions_per_t"])
    return atlas, reads, thetaH, thetaL, lambdas

#%%
def save_mixture(data_file, reads):
    np.save(data_file, reads, allow_pickle=True)

def write_celfie_output(output_file, atlas, atlas_coverage=1000): #no summing
    y = np.hstack([(x*atlas_coverage) for x in atlas])
    y_depths = np.ones((y.shape))
    y_depths.fill(atlas_coverage)
    np.save(output_file, [y, y_depths], allow_pickle=True)

def write_celfie_plus_output(output_file, atlas):
    np.save(output_file, atlas, allow_pickle=True)

def write_epistate_output(output_file, thetaH, thetaL, lambdas):
    np.save(output_file, [thetaH, thetaL, lambdas], allow_pickle=True)

def main(config):
    atlas, reads, thetaH, thetaL, lambdas = generate_data(config)
    save_mixture(config["data_file"], reads)
    for model, outfile in zip(config["models"], config["metadata_files"]):
        if model == "celfie":
            write_celfie_output(outfile, atlas)
        elif model == "celfie-plus":
            write_celfie_plus_output(outfile, atlas)
        else:
            write_epistate_output(outfile, thetaH, thetaL, lambdas)


config = {'m_per_region': 5, 'regions_per_t': 3, 't': 2, "true_alpha": np.array([0.9,0.1]), "coverage": 10,
           "models" :["celfie", "celfie-plus", "epistate"], "data_file":"output.npy",
           "metadata_files":["1.npy", "2.npy", "3.npy"]}
main(config)
#%%