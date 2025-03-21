This folder contains JSON config files that specify various options in cell population simulation configuration.
For each simulation, the user must write an updated config file with options as required and specify the name of the config file using --sim_config argument to the cellpop.py script.
Several config files have been provided as examples.

Default simulation configuration options:

"sim_name": User specified string to identify outputs from a workflow

"time_unit": ["second","minute"], Select the time unit of the single cell model

"exp_time": Integer/float value to specify duration of cell population simulation in hours.

"preinc_time": Integer/float value to specify preincubation duration for individual cells to generate heterogenized cells.

"gen0_time":  Integer/float value to specify simulation time for "gen0" to induce asynchronously cycling starting cell population.

"cc_marker": String identifier for cell cycle marker species.

"cellpop": Integer to specify starting cell population.

"td": Float value to specify anticipated doubling time of the cells in hours,

"mb_tr": Float value to specify maximum possible value of the cell cycle marker trough in original unit.

"mb_peak": Float value to specify minimum height of the cell cycle marker peak,

"tneg_hours": Float value to specify duration in hours from generation 0 appended to generation 1 for successfule detection of cell cycle marker peak,

"timespan_over": Float value to specify duration in hours of additional simulation time to prevent loss of information during downsampling,

"downsample_rate": Integer value to specify downsampling ratio to save disk space, every n-th point is saved output files.

"model_module": Dict to specify options ["load_model","run_model","output"] for single cell model.

{
"load_model": Name of the module that generates required inputs for single cell model. Module and function name must match and must be placed under bin/modules.

"run_model": Name of the module that runs the single cell model. Module and functiona name must match and must be placed under bin/modules,

"output": Array specifying the names of output from "run_model" funciton. "xoutS" and "tout" are essentional, may contain additional elements based on model structure.
}

"stim": Dict to specify growth factor ["egf,"ins","hgf","nrg","pdgf","igf","fgf" ] concentrations in nM.

"apoptosis": Dict to specify options for apoptosis submodel (applicable to SPARCED)

{
      "formula": algebraic formula for apoptosis marker as string.
      "threshold": minimum threshold value for apoptosis marker to define apoptosis event.
}

"drs": Dict to specify options for drug dose response simulations.

{
      "drug": Species name for the drug as string,
      "dose_um": Float value to specify drug dose in uM,
      "rep": String identifier for a replicate,
      "output_partition": Boolean to specify whether outputs for each dose and replicate are saved in individual folder under "sim_name"
}
