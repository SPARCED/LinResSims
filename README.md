***Under construction***


This is a pipeline to generate lineage resolved cell population simulations using the SPARCED single cell model.


Instructions for replicating results from the paper:

1. Verify input files
2. Run scripts/createModel.py to build model
3. Once the model build process is complete, MPI can be used to run cell population simulations using the following command:
    mpirun -np [n_cpu] python cellpop_drs.py --arguments [argument_value]
    Here, n_cpu  is the number of MPI processes that the user decides to use for parallelization.
4. Simulation code accepts the following command line arguments:
    sim_name: An arbitrary string defined by the user to create a directory under sparced/output where simulation outputs will be saved.
    cellpop: An integer specifying the number of starting cells for simulation
    exp_time: Duration of experiment in hours
    drug: String specifying species name for the drug of interest
    dose: Applied concentration of the drug in μM
    egf: Serum EGF concentration in nM
    ins: Serum INS concentration in nM
    hgf: Serum HGF concentration in nM
    nrg: Serum Heregulin concentration in nM
    pdgf: Serum PDGF concentration in nM
    igf: Serum IGF concentration in nM
    fgf: Serum FGF concentration in nM
    Upon completion of simulations, the results are saved to disk in a folder structure corresponding to drug name, replicate identifier and
    drug dose respectively. For a single simulation with a specific replicate of a drug dose, outputs from all cells in each generation are 
    saved in a python pickle object.

5. To generate all drug dose response simulation data, run cellpop_drs.py for:
    * 4 drugs
    * 10 dose levels (including control)
    * 10 replicates of each dose
4. Run analysis scripts to generate inputs to gr_score pipeline
5. Run gr_score python scripts
6. Run plotting script(s)





