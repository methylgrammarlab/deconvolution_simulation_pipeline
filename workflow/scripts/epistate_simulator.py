
import numpy as np
import re


class ReadSimulator:
    '''
    simulate single region with t cell types and m cpgs
    '''
    def __init__(self, t, m, lambda_t, theta_high, theta_low):
        self.t = t
        self.m = m
        self.thetaH = theta_high
        self.thetaL = theta_low
        self.Lt = lambda_t

    def simulate_mixture(self, alpha, coverage):
        mixture = np.ones((coverage,self.m))
        mixture.fill(np.nan)
        #set indicators
        cell_identity = np.random.choice(np.arange(self.t), size=coverage, p=alpha, replace=True)
        mu_identity = []
        for i in range(coverage):
            mu_i = np.random.binomial(1,self.Lt[cell_identity[i]])
            mu_identity.append(mu_i)
            theta = (self.thetaL, self.thetaH)[mu_i]
            mixture[i,:] = np.random.binomial(n=1, p=theta, size=self.m)
        return mixture, cell_identity, np.array(mu_identity)


class BigSimulator:
    '''
    simulate multiple regions
    '''
    def __init__(self, t=25, depth=10, m_per_window=10, windows_per_t=100,
                 thetaL=0.2, thetaH=0.8, lambda_high=1, lambda_low=0, alpha=None):
        self.t = t
        self.depth = depth
        self.m_per_window = m_per_window
        self.windows_per_t = windows_per_t
        self.thetaL = thetaL
        self.thetaH = thetaH
        self.Lhigh = lambda_high
        self.Llow = lambda_low
        self.alpha = alpha

    def get_alpha(self):
        return self.alpha

    def simulate_atlas(self, depth):
        reads = []
        lambdas = []
        z = []
        mu = []

        for cell_type in range(self.t):
            lambda_t = np.ones((self.t))
            lambda_t.fill(self.Lhigh)
            lambda_t[cell_type] = self.Llow #thetaL is indicative of the cell type
            window_maker = ReadSimulator(self.t, self.m_per_window, lambda_t, self.thetaH, self.thetaL)

            for i in range(self.windows_per_t):
                window =[]
                true_z = []
                true_mu = []
                for j in range(self.t):
                    alpha = np.zeros(self.t)
                    alpha.fill(1 / self.t)
                    window_j, true_z_j, true_mu_j = window_maker.simulate_mixture(alpha, depth)
                    window.extend(window_j)
                    true_z.extend(true_z_j)
                    true_mu.extend(true_mu_j)

                reads.append(np.vstack(window))
                lambdas.append(np.hstack(lambda_t))
                z.append(np.hstack(true_z))
                mu.append(np.hstack(true_mu))
        return reads, lambdas, z, mu

    def sample_atlas(atlas, depth):
        '''
        generate reads of pure atlas
        :param atlas: beta values per cell type
        :param depth: number of reads to generate per cell type
        :return: list of reads per atlas region, origin (z values) for reads
        '''
        sample = []
        z = []
        t = atlas[0].shape[0]  # number of cell types
        for i in range(len(atlas)):
            read_source = np.repeat(np.arange(t), depth)
            z.append(read_source)
            probability = atlas[i][read_source, :]
            mixture = np.random.binomial(1, probability)
            sample.append(mixture)
        return sample, z

    def simulate_reads(self):
        reads = []
        lambdas = []
        z = []
        mu = []
        for cell_type in range(self.t):
            lambda_t = np.ones((self.t))
            lambda_t.fill(self.Lhigh)
            lambda_t[cell_type] = self.Llow #thetaL is indicative of the cell type
            window_maker = ReadSimulator(self.t, self.m_per_window, lambda_t, self.thetaH, self.thetaL)
            for i in range(self.windows_per_t):
                window, true_z, true_mu = window_maker.simulate_mixture(self.alpha, self.depth)
                reads.append(window)
                lambdas.append(lambda_t)
                z.append(true_z)
                mu.append(true_mu)
        return reads, lambdas, z, mu

def generate_data(args):
    simulator = BigSimulator(args.t, args.depth, args.m_per_window, args.windows_per_t, args.thetaL, args.thetaH,
                             args.lambda_high, args.lambda_low, args.alpha)
    reads, lambdas, z, mu = simulator.simulate_reads()
    true_alpha = simulator.get_alpha()
    return true_alpha, reads, lambdas, z, mu
#%%
def write_output(reads, lambdas, data_file, metadata_file):
    np.save(data_file, reads, allow_pickle=True)
    np.save(metadata_file, lambdas, allow_pickle=True)

def main(args):
    true_alpha, reads, lambdas, z, mu = generate_data(args)
    write_output(reads, lambdas, args.data_file, args.metadata_file)

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("t", type=int)
parser.add_argument("depth", type=int)
parser.add_argument("m_per_window", type=int)
parser.add_argument("windows_per_t", type=int)
parser.add_argument("thetaL", type=float)
parser.add_argument("thetaH", type=float)
parser.add_argument("lambda_high", type=float)
parser.add_argument("lambda_low", type=float)
parser.add_argument("data_file", type=str)
parser.add_argument("metadata_file", type=str)
parser.add_argument('alpha', type=str)
args = parser.parse_args()
args.alpha = np.array([float(x) for x in re.split(',', args.alpha.strip("[]").replace(" ",""))])
args.alpha = args.alpha / args.alpha.sum()
main(args)
#%%
# true_alpha = np.arange(1,26)
# true_alpha = true_alpha/np.sum(true_alpha)
# simulator = BigSimulator(25, 10, 1, 240, 0.01, 0.99,
#                          1, 0, true_alpha)
# reads, lambdas, z, mu = simulator.simulate_reads()