# Bacterial Slimulator

This repository aims at gathering various SLiM scripts to help people simulate bacterial populations.

It contains the main scripts to run the simulations under a given SLiM model. 

## Dependencies

- python3
    - (msprime)[https://msprime.readthedocs.io/en/latest/]
    - scikit-allel
    - pyslim
- (SLiM 3)[https://messerlab.org/slim/]
- (ms)[https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13]

SLiM and ms should be in your path.

## Quickstart

Clone this repo if you want to share your model, or just download it to reuse an existing one.
Chose a model among the one available under Models/. 
A model defines what we are going to simulate (e.g. constant population size).
Defines the parameters for your model, including whether you want to simulate under WF, nonWF or both, and the number of replicates. 
Examples can be found under Params/.

Then, just run:

`$ bactslim.py path_to_modelparam.json`



