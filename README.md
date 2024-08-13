***Under construction***


This is a pipeline to generate lineage resolved cell population simulations using the SPARCED single cell model.


Instructions for replicating results from the paper:

1. Verify input files
    * Species.txt
    * Ratelaws.txt
    * OmicsData.txt
    * GeneReg.txt
    * StoicMat.txt
2. Run scripts/createModel.py to build model
3. Once the model build process is complete, MPI can be used to run cell population simulations using the following command:
    mpirun -np [n_cpu] python cellpop_drs.py --arguments [argument_value]
    Here, n_cpu  is the number of MPI processes that the user decides to use for parallelization.
4. Simulation code accepts the following command line arguments:
    sim_name: An arbitrary string defined by the user to create a directory under sparced/output where simulation outputs will be saved.
    cellpop: An integer specifying the number of starting cells for simulation
    exp_time: Duration of experiment in hours
    drug: String specifying species name for the drug of interest
    rep: String identifier for the current replicate
    dose: Applied concentration of the drug in μM
    egf: Serum EGF concentration in nM
    ins: Serum INS concentration in nM
    hgf: Serum HGF concentration in nM
    nrg: Serum Heregulin concentration in nM
    pdgf: Serum PDGF concentration in nM
    igf: Serum IGF concentration in nM
    fgf: Serum FGF concentration in nM
    Upon completion of simulations, the results are saved to disk in a folder structure corresponding to drug name, replicate identifier and
    drug dose respectively. For a single simulation with a specific replicate of a drug dose, outputs (temporal species trajectories)
    from all cells in each generation are 
    saved in a python pickle object.
5. To generate all drug dose response simulation data, run cellpop_drs.py for:
    * 4 drugs
    * 10 dose levels (including control)
    * 10 replicates of each dose
6. To generate cell population dynamics (number of alive cells over time) from dose response simulation outputs, run analysis_popdyn.py
    For this, results from all drug dose response simulations in step 5 need to be placed in the "output" folder in the main directory.
    Alternatively, outputs may be placed at a secondary locations and the path must be updated in line 68 of analysis_popdyn.py script.
    Outputs for the cell population dynamics will be saved in the "in_silico_drs_summary" folder underthe output directory.
7. To calculate GR score from the cell population dynamics, input files must be prepared for the gr-score calculation pipeline. For this, step 6
    must be completed first. After this, run analysis_grscore.py to generate the gr-score input file, which will be saved as 'drs_grcalc.tsv'
    in the "in_silico_drs_summary" folder.
8. Take the input file generated at step 7 and run the gr-score calculation pipeline. For this:
    8a. Clone the gr-score git repository: https://github.com/datarail/gr_metrics
    8b. Install all dependencies including Python 2.0
    8c. Go to gr_metrics/SRC/python/scripts
    8d. Run python add_gr_column.py [input_path] > [output_path]    
9. To generate dose response comparison (experiment vs. simulation) plots run analysis_grplots.py after gr scores have been calculated. 





