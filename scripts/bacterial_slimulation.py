#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import os
from subprocess import call
import numpy as np
import sys
import pandas as pd
#import matplotlib.pyplot as plt
import argparse
import logging
import time
#sys.path.append("/home/jean/Documents/ML_genetics/dnadna/")
import summary_statistics as ss
from shutil import copyfile
import distutils.spawn
import recapitate as recap
from multiprocessing import cpu_count, Pool

"""
Script to launch a given slim3 script with WF or nonWF model.
In case of WF, it performs first a burnin with ms, loads the msfile, and start the simulation from there.
In case of nonWF, it performs forward simulation without neutral mutations but keeps track of the tree sequence.
The tree sequence is outputed and loaded in python with pyslim.
Pyslim allows then to recapitate the tree, and msprime will add finally the neutral mutation on the recapitated tree.
"""

def slim(script, run_id, out_dir="."):
    """Function to call slim

    Parameters
    ----------
    script : str
        Path to a slim script

    Returns
    -------
    None

    """
    slim_bin = distutils.spawn.find_executable("slim")
    if slim_bin == None:
        raise OSError("Slim executable not found")
    try:
        with open(os.path.join(out_dir, run_id + ".out"), "w") as outf:
            Slim_cmd = [slim_bin,
                        "-t", # print SLiM's total execution time (in user clock time)
                        "-m", # print SLiM's peak memory usage
                        "-d", """"runId='{}'" """.format(os.path.join(out_dir, run_id)),
                        "-d", """"runIdShort='{}'" """.format(run_id),
                        "-d", """"Ne={}" """.format(Ne),
                        "-d", """"Mu={}" """.format(mutation_rate),
                        "-d", """"Rho={}" """.format(recombination_rate),
                        "-d", """"tractlen={}" """.format(int(tractlen)),
                        "-d", """"genomeSize={}" """.format(int(chr_size)),
                        "-d", """"sampleSize={}" """.format(int(sample_size)),
                        "-d", """"N_generations={}" """.format(int(N_generations)),
                        "-d", """"gcBurnin={}" """.format(gcBurnin),
                        script]
            Slim_cmd2 = " ".join(Slim_cmd)
            logging.debug("Calling Slim:\n{}".format(Slim_cmd2))
            ret = call(Slim_cmd2, stdout=outf, stderr=outf, shell=True)
            if ret != 0:
                logging.error("Slim command: {}\nfailed with error {}".format(Slim_cmd2, ret))
                sys.exit(ret)
    except Exception as e:
        logging.error("Slim failed for run {}: {}".format(run_id, e))
        sys.exit(1)


def run_nonWF(run_id, out_dir_nWF):
    t1 = time.time()
    logging.debug("in run nonWF")

    os.makedirs(out_dir_nWF, exist_ok=True)

    if not os.path.isfile(os.path.join(out_dir_nWF, run_id + ".npz")):
        slim(os.path.join(bactslim_dir,
                          "models",
                          f"{params['model']}",
                          f"nonWF_{params['model']}.slim"),
             run_id,
             out_dir=out_dir_nWF)
        time_nWF_fw = time.time() - t1
        t1 = time.time()
 
        tree = os.path.join(out_dir_nWF, run_id + ".tree")
        snp_mat, pos = recap.recapitate(tree, sample_size, recombination_rate, mutation_rate, Ne, gcBurnin)
        logging.debug("recap done")
        np.savez_compressed(os.path.join(out_dir_nWF, run_id + ".npz"),
                            SNP=snp_mat,
                            POS=pos)
        time_nWF_recap = time.time() - t1

        with open(os.path.join(out_dir_nWF, run_id + ".time"), mode='w') as timefile:
            timefile.write(f"{run_id}\tnonWF\tfwd\t{time_nWF_fw}\n")
            timefile.write(f"{run_id}\tnonWF\tburnin\t{time_nWF_recap}\n")
    return out_dir_nWF

def run_WF(run_id, out_dir_WF):

    os.makedirs(out_dir_WF, exist_ok=True)
    logging.debug("in run WF")
    if not os.path.isfile(os.path.join(out_dir_WF, run_id+".npz")):
        slim(os.path.join(bactslim_dir,
                          "models",
                          f"{params['model']}",
                          f"WF_{params['model']}.slim"),
             run_id,
             out_dir=out_dir_WF)
        ss.convert_ms(os.path.join(out_dir_WF, run_id+".msout"))
        try:
            os.remove(os.path.join(out_dir_WF, run_id+".msout"))
            os.remove(os.path.join(out_dir_WF, run_id+".ms"))
        except FileNotFoundError:
            pass
    else:
        logging.info(f"{run_id} already done, continue...")
    return out_dir_WF

def runner(param_type, run_id, out_dir_nWF, out_dir_WF):

    if param_type in ["both", "nonWF"]:
        logging.info(f"> Starting nonWF simulations for {run_id}")
        out_dir_nWF = run_nonWF(run_id, out_dir_nWF)
        logging.debug(f"< nonWF {run_id} simulations done")

    if param_type in ["both", "WF"]:
        logging.info(f"> Starting WF simulations for {run_id}")
        out_dir_WF = run_WF(run_id, out_dir_WF)
        logging.debug(f"< WF {run_id} simulations done")

def worker(param):
    runner(**param)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("config",
                        help="Path to the config file (json)")

    parser.add_argument('--loglevel',
                         default='INFO',
                         action='store',
                         type=str,
                         metavar='INFO',
                         help='Amount of log info: DEBUG, [INFO], WARNING, ERROR, CRITICAL')

    parser.add_argument('--cpu',
                    help="Number of cpu to run the replicates. By default, all available",
                    type=int,
                    default=cpu_count())

    parser.add_argument('--outdir',
                        help="""Output directory.
                        Results will under [outdir]/model/, where
                        model is defined in the config file.
                        Simulations will not be rerun if the corresponding output file already exists.""",
                        type=str,
                        default=".")


    args = parser.parse_args()
    params = json.loads(open(args.config).read())

    logging.basicConfig(stream=sys.stdout,
                        #filename=os.path.join(out_dir, 'output.log'),
                        level=logging.getLevelName(args.loglevel),
                        format="%(asctime)s; %(levelname)s;  %(message)s",
                        datefmt="%d/%m/%Y %H:%M:%S")

    logging.debug(f"{params['rescaling_factor']}")


    logging.debug(f"rescaling parameters by a factor of {params['rescaling_factor']}")
    recombination_rate = params["rho"] * params["rescaling_factor"]
    logging.debug(f"rho: {params['rho']} -> {recombination_rate}")
    mutation_rate = params["mu"] * params["rescaling_factor"]
    logging.debug(f"mu: {params['mu']} -> {mutation_rate}")
    Ne = int(params["Ne"] / params["rescaling_factor"])
    logging.debug(f"Ne: {params['Ne']} -> {Ne}")

    sample_size = int(params["sample_size"])
    tractlen = params["tract_length"]
    chr_size = int(params["chr_size"])

    try:
        N_generations = int(params["N_generations"] / params["rescaling_factor"])
    except KeyError:
        logging.warning("Number fo generations (N_generations) not set, using 1000")
        N_generations = 1000
    logging.debug(f"Ne: {params['N_generations']} -> {N_generations}")
    if recombination_rate * chr_size > 1:
        logging.warning("""recombination rate too high given chromosome size.
                        >>> *Recombination rate* is modified to 1/chr_size
                        (every indiv will endure gene conversion at each generation)""")
        recombination_rate = 1 / chr_size

    # TODO: Document GC_in_burnin
    try:
        if params["GC_in_burnin"] == 1: # if 1, use actual rec rate
            gcBurnin = recombination_rate # Already rescaled
        elif 0 < params["GC_in_burnin"] < 1: # other rate, rescaled then
            gcBurnin = params["GC_in_burnin"] * params["rescaling_factor"]
        elif params["GC_in_burnin"] > 1:
            logging.error("Param GC_in_burnin should be within [0, 1]")
        else: # if 0 
            gcBurnin = 0
    except KeyError: # not specified
        gcBurnin = 0

    ## Main path
    script_dir = os.path.dirname(os.path.realpath(__file__))
    bactslim_dir = os.path.abspath(os.path.join(script_dir, ".."))
    root_out_dir = os.path.join(args.outdir, params["model"])
    
    ##################
    # Output params  #
    ##################
    os.makedirs(root_out_dir, exist_ok=True)
    copyfile(args.config, os.path.join(root_out_dir,
                                       f"params_{params['model']}-{params['rescaling_factor']}_{params['scenario']}"))



    ##################
    # Output Results #
    ##################
    format_rep = "0>" + str(np.floor(np.log10(params["N_replicat"]) + 1).astype(int))
    df_time = pd.DataFrame(columns=["WF", "nonWF"])
    scen_id = "{model}-{rf}_{scenario}".format(model=params["model"],
                                                         rf=params["rescaling_factor"],
                                                         scenario=params["scenario"])
    out_dir_WF = os.path.join(root_out_dir,
                              "{type}/{scen_id}".format(model=params["model"],
                                                                type="WF",
                                                                scen_id=scen_id))
    out_dir_nWF = os.path.join(root_out_dir,
                                "{type}/{scen_id}".format(model=params["model"],
                                                                  type="nonWF",
                                                                  scen_id=scen_id))

    all_rep_param = []
    for rep in range(1, params["N_replicat"]+1):

        run_id = "{scen_id}_{replicat:{format_rep}}".format(scen_id=scen_id,
                                                             replicat=rep,
                                                             format_rep=format_rep)

        rep_param = {"param_type":params["type"], "run_id":run_id, 
                     "out_dir_nWF":out_dir_nWF, "out_dir_WF":out_dir_WF}
        all_rep_param.append(rep_param)

    logging.info(f"start simulations with {cpu_count()} cpu for {len(all_rep_param)} replicates")

    with Pool(processes=args.cpu) as pool:
        res = pool.map_async(worker, all_rep_param)
        pool.close()
        pool.join()
    logging.info("Done\nComputing sumstats...")

    if params["type"] in ["both", "WF"]:
        ss.do_sum_stats(out_dir_WF, ld_kws={"circular":True}, size_chr=chr_size, label="WF")
    if params["type"] in ["both", "nonWF"]:
        ss.do_sum_stats(out_dir_nWF, ld_kws={"circular":True}, size_chr=chr_size, label="nonWF")
