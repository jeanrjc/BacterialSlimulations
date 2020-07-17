#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
try:
    import allel
except ImportError:
    pass
from scipy.spatial.distance import squareform
import itertools as it
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns
import warnings
import sys
try:
    import ghalton
except ImportError:
    pass

warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', RuntimeWarning)
import re
import logging
import glob
import re
import argparse
from multiprocessing import cpu_count, Pool

def harmonic_nm1(n):
    """return the n-1 harmonic number"""

    return sum([1/(i-1) for i in range(2, n+1)])


def convert_ms(msfile, n_samples=None):
    """
    Read ms file, which contains:
        //
        segsites: 13
        positions: 0.01 0.2 ...
        00100
        00010
    write them on disk in compressed numpy format.
    """

    df = None
    positions = None
    pos_ok = False
    with open(msfile, "r") as file_in:
        sample = 0
        for i, line in enumerate(file_in):
            if "segsites:" in line.split():
                num_segsites = int(line.split()[1])
                if num_segsites == 0:
                    return
            elif "positions:" in line.split():
                positions = np.array([float(j) for j in line.split()[1:]])
                df = pd.DataFrame(columns=range(num_segsites), dtype=int)
                pos_ok = True
            elif pos_ok and line!="\n":
                df.loc[sample] = np.array([int(j) for j in line.split()[0]])
                sample += 1
                if sample == n_samples:
                    break
    if "msin" in msfile:
        outfile = os.path.splitext(msfile)[0] + "_ms"
        # remove columns without snp (because not present in the subsample)
        d = df.values
        # get columns where at least one element in the column is different from the first element
        positions = positions[(d!=d[0]).any(axis=0)]
        d = d[:, (d!=d[0]).any(axis=0)].astype(int)
    else:
        outfile = os.path.splitext(msfile)[0]
        d = df.values.astype(int)
    np.savez_compressed(outfile,
                        SNP=d,
                        POS=positions)

def read_ms_compressed(npzfile, key="all"):
    """
    Takes a .npz file and return all data (SNP and position).
    If one want to get only SNP matrix, set key="SNP",
    or key="POS" for only Position arrays.
    """
    data = np.load(npzfile)
    data = dict(zip((k for k in data), (data[k] for k in data)))
    if any(np.diff(data["POS"]) < 0):
        data["POS"] = np.cumsum(data["POS"])
    if key=="all":
        return data["SNP"], data["POS"]
    else:
        return data[key]

def split_simid(simid):
    """
    Given a sim_ID looking like: "model-A_N_X.extension",
    return the model name (model-A), the scenario number (N) and the replicate number (X)

    Parameters
    ----------
    simid : str
        simulation ID of the form "model-A_N_X.npz".

    Returns
    -------
    model : str
        name of the model.
    scenario : str
        scenario id.
    replicat : int
        replicat number.
    run_id : str
        run id is the N_X part

    Raise
    -----
    NameError :
        If the simulation ID can't be parsed.
    """

    try:
        model, scenario, replicat = re.findall("^([\w\d\-]+)_([\w\d\-]+)_(\d+).*", simid)[0]
    except IndexError:
        raise NameError("The name of the simulation should be of the form `model-A_N_X.extension`")
    return model, scenario, int(replicat), "{}_{}".format(scenario, replicat)

def sfs(haplotype, ac, nindiv=None, folded=False):
    """
    Compute sfs for SNP matrix
    """
    if nindiv == None:
        nindiv = haplotype.shape[1]
    tmp_df = pd.DataFrame({"N_indiv":range(1, nindiv)})
    if folded:
        df_sfs = pd.DataFrame(allel.sfs_folded(ac), columns=["count_SNP"])
        df_sfs["i_xi"] = allel.sfs_folded_scaled(ac)
        df_sfs.index.name = "N_indiv"
        df_sfs.reset_index(inplace=True)
        df_sfs = df_sfs.merge(tmp_df, on="N_indiv", how="right").fillna(0).astype(int)
    else:
        df_sfs = pd.DataFrame(allel.sfs(ac.T[1]), columns=["count_SNP"])
        df_sfs["i_xi"] = allel.sfs_scaled(ac.T[1])
        df_sfs.index.name = "N_indiv"
        df_sfs.reset_index(inplace=True)
        df_sfs = df_sfs.merge(tmp_df, on="N_indiv", how="right").fillna(0).astype(int)

    df_sfs["freq_indiv"] = df_sfs.N_indiv / nindiv
    return df_sfs

def LD(haplotype, pos_vec, size_chr, circular=True, distance_bins=None, gaps_type="short", min_SNP_pairs=300):
    """
    Compute LD for a subset of SNPs drawn with different gap sizes in between them.
    Gap sizes follow power 2 distribution.
    The LD is then computed and averaged over different bin (distance_bins) sizes.

    Parameters
    ----------
    haplotype : numpy 2D array or allel.haplotype
        SNP matrix where in the first dimension are the SNP (rows) and
        in the second dimension (columns) are the samples.
    pos_vec : 1D array
        array of absolute positions in [0, size_chr].
    size_chr : int
        Size of the chromosome.
    circular : bool
        Whether to consider the chromosome circular or not.
        If circular, the maximum distance between 2 SNPs is thus half the chromosome.
    distance_bins : int or list
        LD will be averaged by bins of distances
        e.g. if distance_bins = [0, 100, 1000, 10000], LD will be averaged for the groups [0,100[, [100, 1000[, and [1000, 10000[
        If distance_bins is an int, it defines the number of bins of distances for which to compute the LD
            The bins are created in a logspace
        If distance_bins is a list, they will be used instead
    gaps_type: str
        Pairs of SNP considered are separated by a given number (gap) of columns. Not all pairs are considered.
        By defaut (`short`), gaps are power of 2 up to the closest power of 2 of the number of SNP.
        Meaning that most of the comparisons will be done on close SNPs (short distance).
        If one wants to sample more at large distance (to test for circularity for instance), use `long` instead of `short`
        Using `long` will add gaps like: n_SNP - gaps. It will take more time to run.
    min_SNP_pairs: int
        Minimum number of pairs of SNP to consider for a given gap size.
        If the gap size is big enough such that there is less than `min_SNP_pairs` possible pairs,
        then all pairs are considered.

    Returns
    -------
    DataFrame
        Table with the distance_bins as index, and the mean value of
    """

    if isinstance(distance_bins, type(None)) or isinstance(distance_bins, int):
        if isinstance(distance_bins, int):
            n_bins = distance_bins - 1
        else:
            n_bins = 19
        if circular:
            distance_bins = np.logspace(2, np.log10(size_chr//2), n_bins)
            distance_bins = np.insert(distance_bins, 0, [0])
        else:
            distance_bins = np.logspace(2, np.log10(size_chr), n_bins)
            distance_bins = np.insert(distance_bins, 0, [0])

    n_SNP, n_samples = haplotype.shape

    # gaps are distance between SNPs in term of position in the snp matrix (not in bp)
    gaps_interval = (2 ** np.arange(0, np.log2(n_SNP), 1)).astype(int) # log2 scales of intervals
    if gaps_type.lower() == "long":
        gaps_interval = np.unique(np.concatenate([gaps_interval,
                                                np.array(list(n_SNP//2 - gaps_interval[:len(gaps_interval)//2])[::-1]).astype(int),
                                                np.array(list(n_SNP - gaps_interval)[::-1])])).astype(int)
    else:
        if gaps_type.lower() != "short":
            logging.warning("gaps should be either `short` or `long`. Using short instead of f{gaps_type}")

    selected_snps = []
    for gi, gap in enumerate(gaps_interval):

        if circular:
            max_value = n_SNP
        else:
            max_value = n_SNP - gap
        if max_value < min_SNP_pairs: # min_SNP_pairs : min number of SNP pairs to consider.
            # if not many possible pairs possible, just take them all directly,
            # instead of reaching that number after many more random trials
            snps = np.arange(0, n_SNP, gap)
            snp_pairs = np.unique([((snps[i] + i) % n_SNP, (snps[i + 1] + i) % n_SNP) for i in range(len(snps) - 1)], axis=0)
            snp_pairs = np.concatenate([(snp_pairs + i)%n_SNP  for i in range(max_value)], axis=0)
        else:
            if not circular:
                snps = np.arange(0, n_SNP, gap) + np.random.randint(0, (n_SNP - 1) % gap + 1)  # adding a random start (+1, bc 2nd bound in randint is exlusive)
                # non overlapping contiguous pairs
                # snps=[ 196, 1220, 2244] becomes
                # snp_pairs=[(196, 1220), (1221, 2245)]
                snp_pairs = np.unique([((snps[i] + i) % n_SNP, (snps[i + 1] + i) % n_SNP) for i in range(len(snps) - 1)], axis=0)

                # If we don't have enough pairs (typically when gap is large), we add a random rotation until we have at least 300)
                #count = 0
                # remove pairs that are over the edges
                snp_pairs = snp_pairs[snp_pairs[:, 0] < snp_pairs[:, 1]]
            else:
                snps = np.arange(0, n_SNP, gap) + np.random.randint(0, (n_SNP - 1))  # adding a random start
                # non overlapping contiguous pairs
                # snps=[ 196, 1220, 2244] becomes
                # snp_pairs=[(196, 1220), (1221, 2245)]
                snp_pairs = np.unique([((snps[i] + i) % n_SNP, (snps[i + 1] + i) % n_SNP) for i in range(len(snps) - 1)], axis=0)

            last_pair = snp_pairs[-1]

            while len(snp_pairs) < min(min_SNP_pairs, max_value):
                #count += 1
                #if count % 10 == 0:
                    #print(">>  " + str(gap) + " - " + str(len(np.unique(snp_pairs, axis=0))) + " -- "+ str(len(snps) - 1) + "#" + str(count))
                #remainder = (n_SNP - 1) % gap if (n_SNP - 1) % gap != 0 else (n_SNP - 1) // gap
                shift =  np.random.randint(1, n_SNP) % n_SNP
                new_pair = (last_pair + shift) % n_SNP
                snp_pairs = np.unique(np.concatenate([snp_pairs,
                                                    new_pair.reshape(1, 2) ]), axis=0)
                last_pair = new_pair

                if not circular:
                    snp_pairs = snp_pairs[snp_pairs[:, 0] < snp_pairs[:, 1]]

        selected_snps.append(snp_pairs)

    ld = pd.DataFrame()
    for i, snps_pos in enumerate(selected_snps):

        if circular :
            pos_i = pos_vec[snps_pos]
            min_dist = np.array([min(np.diff(pi)%size_chr, np.diff(pi[::-1])%size_chr) for pi in pos_i])%size_chr/2
            sd = pd.DataFrame(min_dist, columns=["snp_dist"]) # %size_chr/2 because max distance btw 2 SNP is size_chr/2
        else:
            sd = pd.DataFrame((np.diff(pos_vec[snps_pos])), columns=["snp_dist"])

        sd["dist_group"] = pd.cut(sd.snp_dist, bins=distance_bins)
        sr = [allel.rogers_huff_r(snps) ** 2 for snps in haplotype[snps_pos]]
        sd["r2"] = sr
        sd["gap_id"] = i
        ld = pd.concat([ld, sd])

    ld2 = ld.dropna().groupby("dist_group").agg(        
            mean_dist=('snp_dist', 'mean'),        
            mean_r2=('r2','mean'),        
            Count=('r2','count'),        
            sem_r2=('r2','sem') )

    return ld2


def tajimasD(haplotype, pos_vec=None, window=None):
    """
    Given a snp_mat, return tajima's D
    window: Number of window of equal size to slice the position vector.

    if windowed stat, provide pos_vec too.
    """
    if not window:
        allel_count = haplotype.count_alleles()
        return allel.tajima_d(allel_count)
    else:
        all_tajD = []
        pos = pd.DataFrame(pos_vec, columns=["pos"])
        pos["pos_cat"] = pd.cut(pos.pos, 100, labels=range(1, window + 1))
        pos.index.name = "SNP"
        snp_per = pos.reset_index().groupby("pos_cat").SNP.unique()
        for per in snp_per:
            if len(per):
                allel_count = haplotype[per].count_alleles()
                all_tajD.append(allel.tajima_d(allel_count))
            else:
                all_tajD.append(np.nan)
        return pd.DataFrame(all_tajD, columns=["TajD"])


def ihs(haplotype, pos_vec, window=None):
    """Compute the standardize integrated haplotype score"""

    ihs = allel.ihs(haplotype, pos_vec, min_maf=0.01,
                    include_edges=True)
    ihs_stand, bins = allel.standardize_by_allele_count(ihs,
                                                        haplotype.count_alleles().T[1],
                                                        diagnostics=False)
    if window:
            di = pd.DataFrame(ihs_stand, columns=["iHS"])
            di["pos_cat"] = pd.cut(pos_vec, window, labels=range(1, window + 1))
            dig = di.groupby("pos_cat").iHS.mean()
            return dig
    else:
        return ihs_stand

def nsl(haplotype, pos_vec=None, window=None):
    """
    Compute the standardize number of segregating sites by length (nSl)
    for each variant, comparing the reference and alternate alleles,
    after Ferrer-Admetlla et al. (2014)

    if windowed stat, provide pos_vec too.
    """

    nsl = allel.nsl(haplotype)
    nsl_stand, bins = allel.standardize_by_allele_count(nsl,
                                                        haplotype.count_alleles().T[1],
                                                        diagnostics=False)
    if window:
        dn = pd.DataFrame(nsl_stand, columns=["nSL"])
        dn["pos_cat"] = pd.cut(pos_vec, window, labels=range(1, window + 1))
        dng = dn.groupby("pos_cat").nSL.mean()
        return dng
    else:
        return nsl_stand

def worker_do_sum_stats(param):
    do_sum_stats(**param)

def do_sum_stats(scenario_dir, size_chr=2e6,
                 ld_kws=None, sfs_kws=None,
                 label="", nrep="all", overwrite=False):
    """Compute sfs and LD for a set of replicates.

    Parameters
    ----------
    scenario_dir : str
        path to the directory where the outputs (npz file) of the replicates for a given scenario are.
    name_id : str
        identifier for the scenario.
    size_chr : int
        Size of the chromosome.
    ld_kws : dict
        Keywords arguments to pass to the ld function.
        Available kws are: circular[True], distance_bins.
    sfs_kws : dict
        Keywords arguments to pass to the sfs function.
        Available kws are: folded[False].
    label : str
        Give a label for the scenario.
    nrep : int
        Whether to use all replicates (default) or just a subset (for testing purpose).
    overwrite : bool
        If False (default), the output is appended at the end of the existing file.
        Otherwise, it overwrites over it.

    Returns
    -------
    None
        Nothing is returned. Files are written on disk:
        scen_id.mut1
        scen_id.sfs
        scen_id.ld
        scen_id.sel (contains data for Tajima's D, ihs and nsl)
    """

    if ld_kws == None:
        ld_kws = {}
    ld_kws.update({"size_chr":size_chr})
    if sfs_kws == None:
        sfs_kws = {}
    if scenario_dir[-1] == "/":
        scenario_dir = scenario_dir[:-1]
    outdir, scen_id = os.path.split(scenario_dir)
    outdir = os.path.join(outdir, "sumstats")
    os.makedirs(outdir, exist_ok=True)

    all_sfs = pd.DataFrame()
    all_ld = pd.DataFrame()
    all_sel = pd.DataFrame()
    npzfiles = [i for i in os.listdir(scenario_dir) if i.endswith("npz")]
    if nrep != "all":
        npzfiles = npzfiles[:nrep]
    for npzfile in npzfiles:
            snp_mat, pos_vec = read_ms_compressed(os.path.join(scenario_dir, npzfile))
            # convert in total size
            if pos_vec.max() <= 1:
                pos_vec = (pos_vec * size_chr).round().astype(int)
            n_indiv = snp_mat.shape[0]
            haplotype = allel.HaplotypeArray(snp_mat.T)
            allel_count = haplotype.count_alleles()
            derived_allel_count = allel_count.T[1]
            sim_id = os.path.splitext(os.path.basename(npzfile))[0]
            model, scenario, replicat, run_id = split_simid(sim_id)

            try:
                df_sfs = sfs(haplotype, allel_count, **sfs_kws)
                all_sfs = pd.concat([all_sfs, df_sfs])
                # df_sfs.to_csv(os.path.join(scen_id + ".sfs"), sep="\t", index=False, mode="a", header=False)
            except Exception as e:
                logging.error("While computing SFS for {}\n>>> Error: {}".format(sim_id, e))

            try:
                ld = LD(haplotype, pos_vec, **ld_kws)
                ld["sim_id"] = sim_id
                ld["scenario"] = scenario
                ld["run_id"] = run_id
                ld["label"] = label
                all_ld = pd.concat([all_ld, ld])
            except Exception as e:
                logging.error("While computing LD for {}\n>>> Error: {}".format(sim_id, e))

            window = 100
            try:
                taj = tajimasD(haplotype, pos_vec, window=window)
            except Exception as e:
                logging.error("While doing Tajimas'D for {}\n>>> Error: {}".format(sim_id, e))
                taj = pd.DataFrame()
            try:
                ihs_ser = ihs(haplotype, pos_vec, window=window)
            except Exception as e:
                logging.error("While doing IHS for {}\n>>> Error: {}".format(sim_id, e))
                ihs_ser = pd.Series()
            try:
                nsl_ser = nsl(haplotype, pos_vec, window=window)
            except Exception as e:
                logging.error("While doing NSL for {}\n>>> Error: {}".format(sim_id, e))
                nsl_ser = pd.Series()
            try:
                df_sel = pd.concat([taj, ihs_ser, nsl_ser], axis=1)
                df_sel.index.name = "position_percent"
                df_sel.reset_index(inplace=True)
                all_sel = pd.concat([all_sel, df_sel])
            except Exception as e:
                logging.error("While doing concat of selection sumstats for {}\n>>> Error: {}".format(sim_id, e))

            #df_sel.to_csv(os.path.join(scen_id + ".sel"), sep="\t", index=False, mode="a", na_rep="NaN", header=False)

    all_sfs2 = all_sfs.groupby("N_indiv").mean()
    all_sfs2["i_xi_norm"] = all_sfs2.i_xi / all_sfs2.i_xi.mean()
    all_sfs2["freq_indiv"] = all_sfs2.index / n_indiv
    all_sfs2["i_xi_sem_norm"] = all_sfs.groupby("N_indiv").i_xi.sem() / all_sfs2.i_xi.mean()
    all_sfs2["i_xi_std_norm"] = all_sfs.groupby("N_indiv").i_xi.std() / all_sfs2.i_xi.mean()
    all_sfs2["model"] = model
    all_sfs2["scenario"] = scenario
    all_sfs2["label"] = label
    all_sfs2.reset_index(inplace=True)
    writing_mode = "w" if overwrite else "a"
    with open(os.path.join(outdir, scen_id + ".sfs"), writing_mode) as sfsfile:
        all_sfs2.to_csv(sfsfile,
                        sep="\t",
                        index=False,
                        header=False if (sfsfile.tell() and not overwrite) else True)

    all_ld2 = all_ld.groupby("dist_group").mean()
    all_ld2["model"] = model
    all_ld2["scenario"] = scenario
    all_ld2["label"] = label
    all_ld2.reset_index(inplace=True)
    with open(os.path.join(outdir, scen_id + ".ld"), writing_mode) as ldfile:
        all_ld2.to_csv(ldfile,
                       sep="\t",
                       index=False,
                       header=False if (ldfile.tell() and not overwrite) else True)

    all_sel2 = all_sel.groupby("position_percent").mean()
    for c in all_sel2.columns:
        all_sel2[f"{c}_std"] = all_sel.groupby("position_percent")[c].std()
        all_sel2[f"{c}_sem"] = all_sel.groupby("position_percent")[c].sem()
    all_sel2["model"] = model
    all_sel2["scenario"] = scenario
    all_sel2["label"] = label
    all_sel2.reset_index(inplace=True)
    with open(os.path.join(outdir, scen_id + ".sel"), writing_mode) as selfile:
        all_sel2.to_csv(selfile,
                       sep="\t",
                       index=False,
                       header=False if (selfile.tell() and not overwrite) else True)


def load_sum_stats(scen_id, label="", path=""):
    """Load data from scen_id/scen_id.{sfs|ld|sel} and return the 3 df"""

    df_sfs = pd.read_table(os.path.join(outdir, scen_id+ ".sfs"))
    df_ld = pd.read_table(os.path.join(outdir, scen_id+ ".ld"))
    df_sel = pd.read_table(os.path.join(outdir, scen_id+ ".sel"))

    if label != "":
        df_sfs["label"] = label
        df_ld["label"] = label
        df_sel["label"] = label

    return df_sfs, df_ld, df_sel

def plot_fill(x, y, color, ax, step=1, label=None):
    """
    plot the mean as a line and the standard error to the mean as a shade.

    Parameters
    ----------
    x : array
        x-coordinate.
    y : groupby object
        distribution of y-coordinate, from which the mean and sem will be computed.
    color : str
        A color for the mean and the fill area around.
    ax : Axes object
        Where to plot.
    step : int
        to sample the x and y values every `step` values.
    label : str
        Describe the (x, y) curve.

    """
    x = x[::step]
    y_mean, y_sem = y.mean()[::step], y.sem()[::step]
    ax.plot(x, y_mean, color=color, label=label)
    ax.fill_between(x, y_mean + y_sem, y_mean - y_sem, color=color, lw=0, alpha=0.5)


def relative_position(positions, size_chr=0):
    """ Given a numpy array with absolute positions,
        return a numpy array with the relative positions.
        i.e.:
            pos_rel[i] = position[i+1] - position[i]
        if the chromosome is circular:
            pos_rel[0] = (position[-1] - position[0])%size_chr

    Parameters
    ----------
    positions : numpy array
        Absolute positions of SNPs
    size_chr : int
        if 0, the first relative position will be , otherwise it il

    Returns
    -------
    pos_rel: numpy array
        Relative positions of SNPs
    """

    if size_chr == 0:
        return np.ediff1d(positions, to_begin=positions[0])
    else:
        return np.ediff1d(positions, to_begin=(positions[0] - positions[-1]) % size_chr)

def plot_sfs(df_sfs, by="scenario", window=1, ax=None, legend=True):
    """
    Function to plot sfs

    Parameters
    ----------
    df_sfs : dataframe
        Table generated by the sfs() function
    by : str
        Column on which to group the data.
        By default, `scenario`, but it could be `label`, or a another created columns
    step : int
        If you want to subsample your sfs, and use every `step` values to plot
        instead of all.
    ax : matplotlib Axes
        to plot it in a subplot.
    legend: bool
        Whether to plot the legend

    Returns
    -------
    None
    """

    uniq_ID = df_sfs[by].sort_values().unique()
    nsamp = df_sfs.N_indiv.max()
    if len(uniq_ID) > 1:
        dic_color = {j:plt.cm.viridis(int(i*255/(len(uniq_ID)-1))) for i,j in enumerate(uniq_ID)}
    else:
        dic_color = {uniq_ID[0] : plt.cm.viridis(128)}

    if ax == None:
        fig, ax = plt.subplots(1,1)

    for g in df_sfs.groupby(by):
        color = dic_color[g[0]]
        if window > 1:
            g[1].groupby("N_indiv")[["freq_indiv", "i_xi_norm"
             ]].mean().rolling(window).mean().dropna().plot(x="freq_indiv",
                                            y="i_xi_norm",
                                            yerr=g[1].groupby("N_indiv")[["freq_indiv", "i_xi_sem_norm"
                                                 ]].mean().rolling(window).i_xi_sem_norm.mean(),
                                            ax = ax,
                                            label=g[0],
                                            color=color)

        else:
            # It groupby on N_indiv even though it supposed to have 1 value
            # but it depends on the outer group defined by "by=". If it is
            # something else than "scenario", it can have more than 1 value for
            # a given N_indiv category.
            g[1].groupby("N_indiv")[["freq_indiv", "i_xi_norm"
             ]].mean().plot(x="freq_indiv",
                            y="i_xi_norm",
                            yerr=g[1].groupby("N_indiv")[["freq_indiv",
                                                          "i_xi_sem_norm"
                                                        ]].mean(),
                            ax=ax,
                            color=color,
                            label=g[0])

    ax.axhline(1/(nsamp),
                  color="0.5",
                  zorder=10,
                  linestyle="--")
    if legend:
        ax.legend(loc=6, bbox_to_anchor=(1, 0.5))
    else:
        ax.legend_.set_visible(False)

def plot_ld(df_ld, by="scenario", ax=None, legend=True):
    """Function to plot the Linkage Desiquilibrium.

    Parameters
    ----------
    df_ld : DataFrame
        Dataframe generated by the ld() function.
    by : str
        Column on which to group the data.
        By default, `scenario`, but it could be `label`, or a another created columns
    ax : Axes
        to plot in a subplot.
    legend: bool
        Whether to plot the legend

    Returns
    -------
    None
    """

    uniq_ID = df_ld[by].sort_values().unique()
    if len(uniq_ID) > 1:
        dic_color = {j:plt.cm.viridis(int(i*255/(len(uniq_ID)-1))) for i,j in enumerate(uniq_ID)}
    else:
        dic_color = {uniq_ID[0]: plt.cm.viridis(128)}

    if ax == None:
        fig, ax = plt.subplots(1,1)

    for g in df_ld.groupby(by):
        color = dic_color[g[0]]
        g[1].groupby("dist_group").mean().plot(x="mean_dist",
                                               y="mean_r2",
                                               yerr=g[1].groupby("dist_group").mean_r2.sem(),
                                               kind="scatter",
                                               label=g[0],
                                               ax=ax,
                                               color=color,
                                               legend=True)
    ax.set_xscale("log")
    ax.set_xlim(10,
                10**round(np.log10(df_ld.mean_dist.max())))
    if legend:
        ax.legend(loc=6, bbox_to_anchor=(1, 0.5))
    else:
        ax.legend_.set_visible(False)

def change_range(x, min_val, max_val):
    """x is in [0, 1]. Return x in [min_val, max_val]"""

    nx = x * (max_val - min_val)  + min_val
    return nx

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--save_path", default="", type=str)
    parser.add_argument("--simul_params_path", default=None, type=str)
    args = parser.parse_args()
    simul_params = load_dict_from_json(args.simul_params_path)
    os.chdir(args.save_path + '/' + simul_params['model_name'])
    pattern = re.compile('.*' + simul_params['model_name'] + '_|_[0-9].npz$')
    file_paths = glob.glob('*/*.npz')
    name_ids = np.unique([re.sub(pattern, '', f) for f in file_paths])
    for name_id in name_ids:
        do_sum_stats(name_id, simul_params['model_name'], size_chr=simul_params['segment_length'], circular=False)
    df_sfs, df_ld = load_sum_stats(simul_params['model_name'])
    plot_sfs(df_sfs, simul_params['model_name'])
    plt.savefig(simul_params['model_name'] + "_sfs")
    plot_ld(df_ld, simul_params['model_name'])
    plt.savefig(simul_params['model_name'] + "_ld")
