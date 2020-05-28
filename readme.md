# Bacterial Slimulator

This repository aims at gathering various SLiM scripts to help people simulate bacterial populations.

SLiM implements two simulations frameworks, under the Wright-Fisher model or not.
Depending on the framework, the simulation are handled differently.

It contains the main scripts to run the simulations under a given SLiM model. 

## Dependencies

Just create a conda environment by using the following command once you're in this directory:

    conda env create

It will install the following dependencies for you: 

- python 3.7, with in specifically:
    - [msprime](https://msprime.readthedocs.io/en/latest/)
    - [pyslim](https://pyslim.readthedocs.io/)
    - [scikit-allel](scikit-allel.readthedocs.io/)
- [SLiM 3](https://messerlab.org/slim/)
- [ms](https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13)

SLiM and ms should be in your path.

## Quickstart

Clone this repo if you want to share your model, or just download it to reuse an existing one.
Chose a model among the one available under Models/. 
A model defines what we are going to simulate (e.g. a constant population size).
Define the parameters for your model, including whether you want to simulate under WF, nonWF or both, and the number of replicates. 
Examples can be found under params/.

Then, just run:

`$ ./scripts/bacterial_slimulation.py params/ConstantSize.json`

## TODO

- [ ] Clean the code from all stuff unnecessary (in sumstats)
- [ ] add main script ?

