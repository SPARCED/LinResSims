***Under construction***


This is a pipeline to generate lineage resolved cell population simulations using the SPARCED single cell model.

Requirements:

1. OpenMPI for Ubuntu/ MS-MPI for Microsoft Windows
2. Python 3.9 or higher
    Required python packages:
    * amici = 0.11.18
    * mpi4py = 3.1.3
    * scipy = 1.6.2
    * antimony = 2.12.0.3
    * python-libsbml = 5.19



Instructions for replicating results from the paper:

1. Verify input files
    * Species.txt
    * Ratelaws.txt
    * OmicsData.txt
    * GeneReg.txt
    * StoicMat.txt
   For detailed instruction on the input files visit: https://github.com/birtwistlelab/SPARCED
2. Go to scripts directory and run createModel.py to build model. Verify model sbml file (SPARCED.xml) and AMICI-compiled model (SPARCED folder
    in the main directory).
3. Once the model build process is complete, MPI can be used to run cell population simulations using the following command:   
    mpirun -np [n_cpu] python cellpop_drs.py --arguments [argument_value]   
    Here, n_cpu  is the number of MPI processes that the user decides to use for parallelization.
4. Simulation code accepts the following command line arguments:
    * sim_name: An arbitrary string defined by the user to create a directory under sparced/output where simulation outputs will be saved.
    * cellpop: An integer specifying the number of starting cells for simulation
    * exp_time: Duration of experiment in hours
    * drug: String specifying species name for the drug of interest (alpel_EC, nerat_EC, trame_EC, palbo_EC)
    * rep: String identifier for the current replicate
    * dose: Applied concentration of the drug in μM
    * egf: Serum EGF concentration in nM (optional)
    * ins: Serum INS concentration in nM (optional)
    * hgf: Serum HGF concentration in nM (optional)
    * nrg: Serum Heregulin concentration in nM (optional)
    * pdgf: Serum PDGF concentration in nM (optional)
    * igf: Serum IGF concentration in nM (optional)
    * fgf: Serum FGF concentration in nM (optional)   
    For example, to run the first replicate of cell population simulation with 0.003162 &#956 M dose of trametinib for 72 hours, with 100 starting cells, using
    the name 'in_silico_drs' and 16 CPUs, the following command can be entered:
    mpirun -np 16 python cellpop_drs.py --sim_name in_silico_drs --cellpop 100 --exp_time 72 --drug trame_EC --dose 0.003162 --rep rep1   
    Upon completion of simulations, the results are saved to disk in a folder structure corresponding to drug name, replicate identifier and
    drug dose respectively. For a single simulation with a specific replicate of a drug dose, outputs (temporal species trajectories)
    from all cells in each generation are 
    saved in a python pickle object.
5. To generate all drug dose response simulation data, run cellpop_drs.py for:
    * 4 drugs (alpelisib, neratinib, trametinib, palbociclib)
    * 10 dose levels (including control)
    * 10 replicates of each dose
6. To generate cell population dynamics (number of alive cells over time) from dose response simulation outputs, run analysis_popdyn.py
    For this, results from all drug dose response simulations in step 5 need to be placed in the "output" folder in the main directory.
    Alternatively, outputs may be placed at a secondary locations and the path must be updated in line 68 of analysis_popdyn.py script.
    Outputs for the cell population dynamics will be saved in the "in_silico_drs_summary" folder under the output directory.
7. To calculate GR score from the cell population dynamics, input files must be prepared for the gr-score calculation pipeline. For this, step 6
    must be completed first. After this, run analysis_grscore.py to generate the gr-score input file, which will be saved as 'drs_grcalc.tsv'
    in the "in_silico_drs_summary" folder.
8. Take the input file generated at step 7 and run the gr-score calculation pipeline. For this:
    * 8a. Clone the gr-score git repository: https://github.com/datarail/gr_metrics
    * 8b. Install all dependencies including Python 2.0
    * 8c. Go to gr_metrics/SRC/python/scripts
    * 8d. Run python add_gr_column.py [input_path] > [output_path]    
9. Download all experimental dose response datasets (GR-scores) from here: https://www.synapse.org/#!Synapse:syn18456348/ and place them in
    'in_silico_drs_summary/mcf10a_drs_exp'
10. To generate dose response comparison (experiment vs. simulation) plots, copy the gr scores output file from step 8 to the 
    'in_silico_drs_summary' folder and run analysis_grplots.py after gr scores have been calculated.





