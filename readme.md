# Bacterial Slimulator

This repository aims at gathering various SLiM scripts to help people simulate bacterial populations.

SLiM implements two simulations frameworks, under the Wright-Fisher model or not.
Depending on the framework, the simulation are handled differently.

It contains the main scripts to run the simulations under a given SLiM model. 

## Dependencies

Just create a conda environment by using the following command once you're in this directory:

    conda env create

It will install the following dependencies for you: 

- python 3.7, with, specifically:
    - [msprime](https://msprime.readthedocs.io/en/latest/)
    - [pyslim](https://pyslim.readthedocs.io/)
    - [scikit-allel](scikit-allel.readthedocs.io/)
- [SLiM 3](https://messerlab.org/slim/)
- [ms](https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13)


## Quickstart

Clone this repo if you want to share your model, or just download it to reuse an existing one.
Chose a model among the one available under Models/. 
A model defines what we are going to simulate (e.g. a constant population size).
Define the parameters for your model, including whether you want to simulate under WF, nonWF or both, and the number of replicates. 
Examples can be found under params/.

Then, just run:

    $ ./scripts/bacterial_slimulation.py params/ConstantSize.json

It will create a folder in your current directory named `ConstantSize` (after the `model`'s parameter in the `.json`'s file) with the following architecture :

    ConstantSize
    ├── nonWF
    │   ├── ConstantSize-100_myScenario
    │   │   ├── ConstantSize-100_myScenario_1.npz
    │   │   ├── ConstantSize-100_myScenario_1.out
    │   │   ├── ConstantSize-100_myScenario_1.time
    │   │   ├── ConstantSize-100_myScenario_1_final.trees
    │   │   ├── ...
    │   └── sumstats
    │       ├── ConstantSize-100_myScenario.ld
    │       ├── ConstantSize-100_myScenario.sfs
    └── WF
        ├── ConstantSize-100_myScenario
        │   ├── ConstantSize-100_myScenario_1.npz
        │   ├── ConstantSize-100_myScenario_1.out
        │   ├── ConstantSize-100_myScenario_1.time
        │   ├── ...
        └── sumstats
            ├── ConstantSize-100_myScenario.ld
            ├── ConstantSize-100_myScenario.sfs

With a subfolder for WF and nonWF types of simulation. 
Within each of these folders, there will be one folder per scenario (as per the `scenario:` line in the parameter's file.) 
It is thus possible to simulate different scenario within the same model, and they will be stored in the same model's scenario. 
Replicates are available within each scenario's folder and contains different types of file:
- `.out`: the output of SLiM
- `.npz`: the snp matrix an position vector in numpy compressed file format
- `.time`: the timing of the burnin and forward part
- `_final.trees`: The tree sequence after recapitation, in the nonWF version. 

