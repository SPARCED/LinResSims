# LinResSims User Guide

_Written by Jonah R. Huggins, Arnab Mutsuddy, & Aurore Amrit_

## Overview

This is a pipeline to generate lineage resolved cell population simulations using the SPARCED single cell model.

**LinRessims** runs seamlessly on **Ubuntu 22.04 LTS** , either as a virtual machine (i.e. **VirtualBox**), a container (**Singularity** or **Dockerfile**), or on **Native** **Linux.** This guide should work even if you are using another hypervisor than VirtualBox or that you are running Ubuntu directly on your computer. With a few arrangements, the described steps should also work for other versions of Ubuntu or any Debian-based Linux distribution.

## Installation Guide

For most users, we highly recommend using one of the provided container options. This provides an OS-agnostic approach for users to interact with LinResSims simulation tools without having to hassle with the SPARCED package dependencies. For HPC use, see Singularity Installation.

### Docker Installation

For users with administrator access, (i.e. Linux, MacOS and Windows), we strongly recommend pulling an image of the Docker container ([install Docker here](*[Docker](https://docs.docker.com/get-started/get-docker/))).

**Preferred installation:**

```
docker pull birtwistlelab/linressims:latest
```

Alternatively, the Docker container can be built locally in the event changes are made to the source code:

```
docker buildx build -t birtwistlelab/linressims -f /path/to/LinResSims/container/Dockerfile .
```

Congratulations! You now have a full setup of LinResSims! 🦠

### Singularity Installation

[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) is a containerization platform designed specifically for high-performance computing (HPC) and research environments. It allows users to create, distribute, and execute portable, reproducible containers across different systems. Unlike Docker, Singularity focuses on usability in environments where users don't have root access, such as shared HPC clusters.

#### Building the Singularity Container from a Definition File

To build a container using the `linressims.def` file, make the following alterations to the `container/linressims.def` file:

1. On line 6, specify the absolute path of the host system to the LinResSims directory (e.g. ****`/home/username/LinResSims`****)
2. On line 49, specify the version of OpenMPI running on the host system.
   * Note, this build process has only been tested with OpenMPI. Using wget, the definition file pulls a **user-specified** version of OpenMPI and builds from source within the container following the [instructions specified here](https://apptainer.org/docs/user/1.0/mpi.html#open-mpi-hybrid-container). If the host system and container versions of MPI do not match, this code will not work as intended (e.g. `export OMPI_VERSION=5.0.1`)

After the above alterations are made, execute the following command from the project root directory:

```bash
singularity build --fakeroot container/linressims.sif container/linressims.def
```

* **`linressims.sif`:** The output file, a Singularity Image File (SIF).
* **`linressims.def`:** The definition file that specifies the container's environment and setup.
* **`--fakeroot`:** Flag for building the singularity container without root access

To verify that the container is successfully built, execute the following:

```bash
singularity inspect container/linressims.sif
```

Congratulations! You now have a full setup of LinResSims! 🦠

### Ubuntu 22.04 (local) Installation

For users with administrator (root) access who want to install project dependencies locally outside of a container, an installation script has been provided (`LinResSims/install.sh`) to simplify the dependency installations. To run, execute the following commands:

```
chmod +x ./install.sh
./install.sh
```

Congratulations! You now have a full setup of LinResSims! 🦠

---

## Operation

Operating the LinResSims code can be done either within a container, outside of a container, or at the command line locally. **Execution of the LinResSims code must be performed from within the `LinResSims/scripts/` directory.**

### Overview

#### SPARCED Compilation

By default, each single cell in a population is simulated using the SPARCED model. Our use of the ODE solver **[AMICI](https://amici.readthedocs.io/en/v0.11.12/)** necessitates the model be re-constructed in a C++ directory relative to the SBML path. Therefore, the model compilation step must be executed before it can be run in a python environment. 

- Models not using the AMICI simulator are not subject to this constraint (see `bin/modules/RunTyson.py` for an example)

 To compile the SPARCED model, the user must change directory to /scripts and run the following command:

```
python createModel.py
```

* Compilation takes several minutes to run and provides sparse output while executing.
* An SBML file (SPARCED.xml) and AMICI-compiled model (SPARCED folder in the main directory) serve as verification that this step was successful.

⚠️SPARCED compilation is only necessary:

* Before simulating LinResSims for the first time
* When modifying the input files

#### Code Execution

To run simulations, execute the following command:

```
mpirun -n <CORES> python cellpop.py --sim_config <name_of_config_file>
```

**Flags:**

`-n`: An MPI-specific flag for defining how many processor cores are utilized by a particular instance of code. Default value is 1 core.

`--sim_config`: Specifies the name of the simulation configuration file should be used for simulation with the LinRessims code.

### Interactive Execution

#### Docker

To run the LinResSims tool from interactively within the Docker container, execute the following:

```
docker run -it birtwistlelab/linressims:latest
```

Users are further able to bind the local LinResSims directory with the container LinResSims directory, enabling native, local modifications on the host system:

```
docker run -it --rm -v </path/to>/LinResSims:/LinResSims birtwistlelab/linressims:latest
```

**Flags:**

* **`--rm`** (*Remove*): Automatically removes the container when it stops to prevent the accumulation of stopped containers that would otherwise take up system resources.
* **`-i`** (*Interactive*): Keeps the standard input (`stdin`) open, even if not attached to a terminal. Allows the container to accept input from the user during runtime.
* * Especially useful when paired with `-t` for running an interactive shell session.
* **`-t`** (*TTY, teletypewriter*): Allocates a pseudo-terminal for the container. Allows for better interactivity, like being able to run `bash` or `sh` inside the container and see a command prompt. Often paired with `-i` for fully interactive sessions.
* **`-v`** (*Volume*): Binds a directory to the container volume. In simpler terms, it allows users to link a local directory to a container directory, enabling files to be shared between. Allows for seamless operation of LinResSims on a users personal device with full file sharing as if the tool was installed locally.

#### Singularity

To open a shell inside the container with the `LinResSims` directory bound:

```bash
singularity shell --bind /path/to/host/LinResSims:/LinResSims container/linressims.sif
```

* **`/path/to/host/LinResSims`:** Replace this with the absolute path to your host's LinResSims directory.
* **`/LinResSims`:** This is the directory inside the container where the host directory will be accessible.

**Flags:**

* **`--bind`** : Option to link a host directory into the container

Once inside the container, you'll see a prompt. By default the container launches in the LinResSims directory:

```bash
$ pwd
/LinResSims # Output
```

**Exit the Container:**
Type `exit` to leave the container.

### Batch Execution

Executing the LinResSims code as a batch script (which is necessary for most HPC job schedulers) is also possible.

#### Docker

To run the LinResSims code from outside of the Docker container, execute the following:

```
docker run --rm -v <\path\to\>LinResSims:/LinResSims birtwistlelab/linressims:latest bash -c "
cd scripts
mpirun -n <CORES> python cellpop.py --sim_config <name_of_config_file>
""
```

Here, a bash command is passed to the container to change the directory to scripts, then execute the code.

#### Singularity

To run the LinResSims code from outside of the Singularity container, execute the following:

```
mpirun -n <CORES> singularity exec container/linressims.sif bash -c "
cd scripts
python cellpop.py --sim_config <name_of_config_file> 2>&1 | tee cellpop.log
"
```

Again, a bash command is passed to the container on execution to perform both operations within a single command instance. The Additional `2>&1 | tee cellpop.log` redirects information sent to standard output and standard error to the file `cellpop.log`. This is entirely optional to include.

To demonstrate running the singularity container on an HPC system with SLURM job scheduler, examples of batch scripts have been provided for creating the container (`slurm_files/new-container.sh`), compiling an instance of the SPARCED model (`slurm_files/compile-container.sh`), and for simulating a cell population (`slurm_files/run-container.sh`).

### Additional Simulation Flags

To override the configuration file without writing over existing simulation settings, use the following (optional) command line arguments:

* `--sim_name`: An arbitrary string defined by the user to create a directory under sparced/output where simulation outputs will be saved.
* `--cellpop`: An integer specifying the number of starting cells for simulation
* `--exp_time`: Duration of experiment in hours
* `--rep`: String identifier for the current replicate
* `--egf`: Serum EGF concentration in nM
* `--ins`: Serum INS concentration in nM
* `--hgf`: Serum HGF concentration in nM
* `--nrg`: Serum Heregulin concentration in nM
* `--pdgf`: Serum PDGF concentration in nM
* `--igf`: Serum IGF concentration in nM
* `--fgf`: Serum FGF concentration in nM

## Configuration Files

Workflow variables for cell population simulations are specified with the use of a json configuration file, which the user may define for each simulation run. This allows the alteration of several key workflow parameters without modification of the simulation script itself. By default simulation config files are located in the folder `LinResSims/sim_configs/`. For a detailed overview of the structure and keys of the configuration file, see `LinResSims/sim_configs/README.md`

## Reproducing Published Examples

To simplify reproduction of our results,  bash scripts ()executable on a SLURM job scheduler) have been provided at `LinResSims/slurm_files`. Please execute these scripts in the following order:

1. `LinResSims/slurm_files/new-container.sh # Builds the singularity container on the local system.`
2. `LinResSims/slurm_files/compile-container.sh #builds an AMICI model for SPARCED simulation`
3. `LinResSims/slurm_files/run-container.sh # Executes a single simulation of the SPARCED model based on the default_SPARCED.json settings`  OR
   `LinResSims/slurm_files/figure_2defg.sh # Runs the simulations necessary to reproduce Figures 2D-F.`

As a demonstration,   an example command to run the first replicate of cell population simulation with 0.003162 μM dose of trametinib for 72 hours, with 100 starting cells, using
the name 'in_silico_drs' and 16 CPUs is provided below:

```
mpirun -n 16 python cellpop.py --sim_name in_silico_drs --cellpop 100 --exp_time 72 --drug trame_EC --dose 0.003162 --rep rep1
```

Upon completion of simulations, the results are saved to disk in a folder structure corresponding to drug name, replicate identifier and drug dose respectively (e.g. `LinResSims/output/in_silico_drs/drs_trame/drs_trame_rep1/trame_EC_0.003162/` ). For a single simulation with a specific replicate of a drug dose, outputs (temporal species trajectories) from all cells in each generation are saved in a python pickle object (e.g. `LinResSims/output/in_silico_drs/drs_trame/drs_trame_rep1/trame_EC_0.003162/output_g1.pkl`).

### Visualization

To replicate the figures from the paper that use simulation outputs, dose resopnse simulaitons for all 4 drugs, across 10 specified dose levels and 10 replicates must have been completed using a unique "sim_name" ("in_siloco_drs" by default) and placed at a convenient location (LinResSims/output by default). For making plots using simulation outputs from a given drug dose and replicate, we have provided a python class "drs_dict" defined within "bin/modules/drsPlotting.py". Example of its usage to generate several types of plots have been provided as jupyter notebooks under LinResSims/jupyter_notebooks.

* figure_1c.ipynb: cross generational protein level trajectories and single cell lineage tree
* figure_2abc.ipynb: cell population dendrogram with control and dosage populations.

Some population level visualizations rely on cell population dynamics and require further analysis after simulation. For example, cell population dynamics require alive cell counts over time to have been completed. 

To generate cell population dynamics (number of alive cells over time) from dose response simulation outputs, run analysis_popdyn.py. For this, results from all drug dose response simulations need to be placed in the "output" folder in the main directory. Alternatively, outputs may be placed at a secondary locations and the path must be updated in line 68 of analysis_popdyn.py script. Outputs for the cell population dynamics will be saved in the "in_silico_drs_summary" folder under the output directory.

Furthermore, visualizing dose response for mutiple drugs, doses, and replicates in terms of GR-score, requires the calculation of GR score after the cell population dynamics have been computed. Further instructions have been provided in the following section. 

### Calculating GR Scores

To calculate GR score from the cell population dynamics, input files must be prepared for the gr-score calculation pipeline. The below steps describe calculating GR scores from results:

1. Ensure that script `figure_2defg.sh` as been executed
2. Complete the **Visualization** instructions provided in the previous section.
3. Run `analysis_grscore.py` to generate the gr-score input file, which will be saved as `drs_grcalc.tsv` in the `in_silico_drs_summary` folder.
4. Take the input file generated at step 3 and run the gr-score calculation pipeline:
   1. Clone the gr-score git repository:` https://github.com/datarail/gr_metrics`
   2. Install all dependencies including Python 2.0
   3. Go to `gr_metrics/SRC/python/scripts`
   4. Run `python add_gr_column.py [input_path] > [output_path]`
5. Download all experimental dose response datasets (GR-scores) from here: https://www.synapse.org/#!Synapse:syn18456348/ and place them in `in_silico_drs_summary/mcf10a_drs_exp`
6. Plots from figures 1,2 can be generated with jupyter notebooks included in the `LinResSims/jupyter_notebooks` folder.

## Running cell population simulations with a new single cell model:

By default, the cell population simulation workflow uses the SPARCED single cell model. It is capable of running simulations with a different single cell model given that the model has a compatible structure. A compatible model must satisfy the following requirements:

* The model must have a state matrix representing a single cell.
* The model must have a variable representing dynamic molecular signature of cell cycle markers, i.e., periodic activation and inactivation of cyclins.
* The model must be executable within a python module.

To replace the SPARCED model in cell population simulations with another single cell model:

1. Place all single cell simulation operations within a python function (see `LinResSims/bin/modules/RunTyson.py` for an example).
2. Write another python function to generate an input dict for the single cell model function, mirroring the input/output structure of the LoadSPARCED function (see `LinResSims/bin/modules/LoadTyson.py` for an example).
3. Save both python functions as modules with the same name as the functions under `LinResSims/bin/modules`.
4. Write a json config file with key-specific values appropriate for the new model structure. Be sure to make "load_model" and "run_model" options consistent with the new module names. For more details on the stucture of the sim config, see `sim_configs/README.md`

The Tyson 1991 cell cycle model has been presented as an example for this procedure. The "load_model" and "run_model" modules have been provided as `LinResSims/bin/modules/LoadTyson.py `and `LinResSims/bin/modules/RunTyson.py`. The sim_config json file corresponding to this workflow is ` LinResSims/sim_config/default.json`

## Contributors Guide

To enable broader portability of the LinResSims project, source code and dependencies are packaged into a distributable wheel using the `pyproject.toml` file and python's `build` command. Further, packages are installed at `/usr/local/lib/python3.10/site-packages/` via the `pip` package manager. In the event that you wish to contribute to update or change python packages, the distributable files (located at `LinResSims/dist/`) will need to be updated as well for the changes to take affect.

1. Update the `pyproject.toml` file with any modifications to the package lists under the `dependencies` variable (line 18)
2. From the project root directory, execute the following command
   * `python -m build`
3. Install the new packages using pip:
   * `pipinstalldist/linressim-1.0-py3-none-any.whl --verbose --force`
