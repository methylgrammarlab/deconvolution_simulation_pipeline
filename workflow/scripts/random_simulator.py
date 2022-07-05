#################################################
# FILE: random_simulator.py
# WRITER: Irene Unterman
# DESCRIPTION: simulate methylation sequencing
# with an atlas of random methylation probabilities
#
#################################################
import numpy as np
import re

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

def generate_mixture(atlas, alpha, depth):
    '''
    sample WGBS sequences from an atlas of methylation probabilities
    :param atlas: t cell types by m CpGs. filled with beta values
    :param alpha: proportions of cell types to sample. ordered by atlas
    :param depth: amount of reads to sample, each will be the entire
                    length of the atlas (m CpGs)
    :return: depthXm methylayion np array
    '''
    t, m = atlas.shape #cell types, cpgs
    z = np.random.choice(
        t,
        depth,
        replace=True,
        p=alpha,
    )  # assign reads based on the cell type proportions

    probability = atlas[z, :]

    mixture = np.random.binomial(1, probability)

    return mixture

def sample_atlas(atlas, depth):
    '''
    generate reads of pure atlas
    :param atlas: beta values per cell type
    :param depth: number of reads to generate per cell type
    :return: list of reads per atlas region, origin (z values) for reads
    '''
    sample = []
    z = []
    t = atlas[0].shape[0] #number of cell types
    for i in range(len(atlas)):
        read_source = np.repeat(np.arange(t), depth)
        z.append(read_source)
        probability = atlas[i][read_source,:]
        mixture = np.random.binomial(1, probability)
        sample.append(mixture)
    return sample, z

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


def generate_data(args):
    atlas = generate_atlas(args.t, args.m_per_window, args.windows_per_t*args.t)
    reads = [generate_mixture(x, args.alpha, args.depth) for x in atlas]
    thetaH, thetaL, lambdas =  make_thetas_and_lambdas(atlas, args.t, args.windows_per_t)
    return atlas, reads, thetaH, thetaL, lambdas
#%%
def write_output(reads, atlas, data_file, thetaH, thetaL, lambdas, celfie_metadata_file, epistate_metadata_file):
    np.save(data_file, reads, allow_pickle=True)
    np.save(celfie_metadata_file, atlas, allow_pickle=True)
    np.save(epistate_metadata_file, [thetaH, thetaL, lambdas], allow_pickle=True)


def main(args):
    atlas, reads, thetaH, thetaL, lambdas = generate_data(args)
    write_output(reads, atlas, args.data_file, thetaH, thetaL, lambdas,
                 args.celfie_metadata_file, args.epistate_metadata_file)

import argparse


parser = argparse.ArgumentParser()
parser.add_argument("t", type=int)
parser.add_argument("depth", type=int)
parser.add_argument("m_per_window", type=int)
parser.add_argument("windows_per_t", type=int)
parser.add_argument("data_file", type=str)
parser.add_argument("celfie_metadata_file", type=str)
parser.add_argument("epistate_metadata_file", type=str)
parser.add_argument('alpha', type=str)
args = parser.parse_args()
args.alpha = np.array([float(x) for x in re.split(',', args.alpha.strip("[]").replace(" ",""))])
args.alpha = args.alpha / args.alpha.sum()
main(args)
#%%
