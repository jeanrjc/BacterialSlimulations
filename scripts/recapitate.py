#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import os
from subprocess import call
import numpy as np
import sys
import pandas as pd
import msprime
import pyslim
import logging

#import matplotlib.pyplot as plt
"""
Script to recapitate haploid simulation
"""


def sample_treeseq(ts, n):
    """Sample n individuals from the TreeSequence.

    Parameters
    ----------
    ts : SlimTreeSequence
        Full TreeSequence from fwd simulation with tree recording
    n : int
        Size of the sample

    Returns
    -------
    SlimTreeSequence
        The subsample tree sequence
    """

    samples = np.random.choice(ts.samples(),
                               size=n,
                               replace=False) # choose n random leaves
    ts_samp = ts.simplify(samples=samples.astype("int32"))
    ts_samp_sts = pyslim.SlimTreeSequence(ts_samp)
    return ts_samp_sts

def recapitate(treesfile, sample_size, recombination_rate, mutation_rate, Ne):

        # load trees with pyslim
        ts = pyslim.load(treesfile)
        num_individuals_0 = len(ts.individuals_alive_at(0))
        n_roots = pd.Series([t.num_roots for t in ts.trees()]).value_counts().to_frame(name="num_tree_with_num_roots")
        logging.debug(f"""The tree sequence has {ts.num_trees} trees on a genome of length {ts.sequence_length}
                         {num_individuals_0} alive individuals, {ts.num_samples} 'sample' genomes
                         and {ts.num_mutations} mutations.
                         number of roots per tree: 
                         {n_roots.__str__()[:-12]}""")

        # discard second genomes (diploids) from the tree.
        #ts_haploid = ts.simplify(samples=[ind.nodes[0] for ind in ts.individuals()])
        ts_recap = ts.recapitate(recombination_rate=1e-20, 
                                 Ne=Ne)
                                 #population_configurations=[msprime.PopulationConiguration(initial_size=Ne)])
        
        # simplify to a subset of the haploids
        sample_inds = np.random.choice(ts_recap.individuals_alive_at(0),
                                       size=sample_size,
                                       replace=False) # choose n random leaves                                       
        sample_nodes = [ts_recap.individual(i).nodes[0] for i in sample_inds]
        ts_samp = ts_recap.simplify(samples=sample_nodes)
        
        n_roots = pd.Series([t.num_roots for t in ts_samp.trees()]).value_counts().to_frame(name="num_tree_with_num_roots")
        logging.debug(f"""The tree sequence has {ts_samp.num_trees} trees on a genome of length {ts_samp.sequence_length}
                         {ts_samp.num_individuals} alive individuals, {ts_samp.num_samples} 'sample' genomes
                         and {ts_samp.num_mutations} mutations.
                         number of roots per tree: 
                         {n_roots.__str__()[:-12]}""")

        # mutate
        ts_mutated = pyslim.SlimTreeSequence(
                        msprime.mutate(ts_samp,
                                    rate=mutation_rate/2, # To have 2.Ne.mu and not 4.Ne.mu
                                    keep=True) # keep existing mutations
                        )
        genotype_matrix = ts_mutated.genotype_matrix()
        snp_mat = genotype_matrix.T
        pos = np.round(ts_mutated.tables.asdict()["sites"]["position"]).astype(int)

        return snp_mat, pos

# # Without actual recapitation
# def fake_recapitate(treesfile, sample_size, recombination_rate, mutation_rate, Ne):

#         # load trees with pyslim
#         ts = pyslim.load(treesfile)
#         num_individuals_0 = len(ts.individuals_alive_at(0))
#         logging.debug(f"""The tree sequence has {ts.num_trees} trees on a genome of length {ts.sequence_length}
#                          {num_individuals_0} alive individuals, {ts.num_samples} 'sample' genomes
#                          and {ts.num_mutations} mutations.""")

#         # discard second genomes (diploids) from the tree.
#         ts_haploid = ts.simplify(samples=[ind.nodes[0] for ind in ts.individuals()])
#         #ts_recap = ts.recapitate(recombination_rate=1e-20, 
#          #                        Ne=Ne)
#                                  #population_configurations=[msprime.PopulationConiguration(initial_size=Ne)])

#         # Keep only alive individuals
#         ts_alive = ts_haploid.simplify(samples=[ind for ind in ts_haploid.individuals_alive_at(0)])
#         logging.debug(f"""The happloid tree sequence has {ts_alive.num_trees} trees on a genome of length {ts_alive.sequence_length}
#                      {ts_alive.num_individuals} alive individuals, {ts_alive.num_samples} 'sample' genomes
#                      and {ts_alive.num_mutations} mutations.""")

#         # subsample the number of individuals.
#         ts_sample = sample_treeseq(ts_alive, sample_size)

#         ts_mutated = msprime.mutate(ts_sample,
#                                     rate=mutation_rate/2, # To have 2.Ne.mu and not 4.Ne.mu
#                                     keep=True) # keep existing mutations

#         genotype_matrix = ts_mutated.genotype_matrix()
#         snp_mat = genotype_matrix.T
#         pos = np.round(ts_mutated.tables.asdict()["sites"]["position"]).astype(int)

#         return snp_mat, pos


if __name__ == "__main__":
    pass
