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

To run the LinResSims tool from within the container, execute the following:

```
docker run -it birtwistlelab/linressims:latest
```

Further, users are able to bind their local LinResSims directory with the container's LinResSims directory, allowing for native, local modifications on the host system:

```
docker run -it --rm -v </path/to>/LinResSims:/LinResSims birtwistlelab/linressims:latest
```

**Flags:**

* **`--rm`** (*Remove*): Automatically removes the container when it stops to prevent the accumulation of stopped containers that would otherwise take up system resources.
* **`-i`** (*Interactive*): Keeps the standard input (`stdin`) open, even if not attached to a terminal. Allows the container to accept input from the user during runtime.
* * Especially useful when paired with `-t` for running an interactive shell session.
* **`-t`** (*TTY, teletypewriter*): Allocates a pseudo-terminal for the container. Allows for better interactivity, like being able to run `bash` or `sh` inside the container and see a command prompt. Often paired with `-i` for fully interactive sessions.
* **`-v`** (*Volume*): Binds a directory to the container volume. In simpler terms, it allows users to link a local directory to a container directory, enabling files to be shared between. Allows for seamless operation of LinResSims on a users personal device with full file sharing as if the tool was installed locally.

Congratulations! You now have a full setup of LinResSims! 🦠

### Singularity Installation

[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) is a containerization platform designed specifically for high-performance computing (HPC) and research environments. It allows users to create, distribute, and execute portable, reproducible containers across different systems. Unlike Docker, Singularity focuses on usability in environments where users don't have root access, such as shared HPC clusters. 

#### Building the Singularity Container from a Definition File

To build a container using the `sparced.def` file, make the following alterations to the `container/linressims.def` file:

1. On line 6, specify the absolute path of the host system to the LinResSims directory (e.g. ****`/home/username/LinResSims`****)
2. On line 49, specify the version of OpenMPI running on the host system.
   1. Note, this build process has only been tested with OpenMPI. Using wget, the definition file pulls a **user-specified** version of OpenMPI and builds from source within the container following the [instructions specified here](https://apptainer.org/docs/user/1.0/mpi.html#open-mpi-hybrid-container). If the host system and container versions of MPI do not match, this code will not work as intended.

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

#### Working Within the Singularity Container

1. **Bind the Host Directory:**
   Use the `--bind` option to link a host directory into the container.
2. **Run the Container:**
   To open a shell inside the container with the `LinResSims` directory bound:

   ```bash
   singularity shell --bind /path/to/host/LinResSims:/LinResSims container/linressims.sif
   ```

   * **`/path/to/host/LinResSims`:** Replace this with the absolute path to your host's LinResSims directory.
   * **`/LinResSims`:** This is the directory inside the container where the host directory will be accessible.
3. **Operate Within the Container:**
   Once inside the container, you'll see a prompt. Run commands as if you're working on a standalone system:

   ```bash
   python createModel.py
   ```
4. **Exit the Container:**
   Type `exit` to leave the container.

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

### Code Execution Overview

To run instances of the SPARCED model, model compilation must first be performed. This is not necessary for other models (e.g. the Tyson model):

```
python createModel.py
```

* Verify model compilation via the production of the sbml file (SPARCED.xml) and AMICI-compiled model (SPARCED folder in the main directory).

To operate, simply pass the following command, either within a container or at the command line:

```
mpirun -n <CORES> python cellpop.py --sim_config <name_of_config_file>
```

**Flags:**

`-n`: An MPI-specific flag for defining how many processor cores are utilized by a particular instance of code. Default value is 1 core.

`--sim_config`: Specifies which simulation configuration file should be used for simulation with the LinRessims code. More on this in the **Configuration File** subsection.

### Container Execution

Executing the LinResSims code from outside the container (which is necessary for most HPC job schedulers) is also possible.

#### Docker Execution

To run the LinResSims code from outside of the Docker container, execute the following:

```
docker run --rm -v <\path\to\>LinResSims:/LinResSims birtwistlelab/linressims:latest bash -c "
cd scripts
mpirun -n <CORES> python cellpop.py --sim_config <name_of_config_file>
""
```

Here, a bash command is passed to the container to change the directory to scripts, then execute the code.

#### Singularity Execution

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

Simulation code accepts the following (optional) command line arguments:

* `--sim_name`: An arbitrary string defined by the user to create a directory under sparced/output where simulation outputs will be saved.
* `--cellpop`: An integer specifying the number of starting cells for simulation
* `--exp_time`: Duration of experiment in hours
* `--drug`: String specifying species name for the drug of interest (alpel_EC, nerat_EC, trame_EC, palbo_EC)
* `--rep`: String identifier for the current replicate
* `--dose`: Applied concentration of the drug in μM
* `--egf`: Serum EGF concentration in nM
* `--ins`: Serum INS concentration in nM
* `--hgf`: Serum HGF concentration in nM
* `--nrg`: Serum Heregulin concentration in nM
* `--pdgf`: Serum PDGF concentration in nM
* `--igf`: Serum IGF concentration in nM
* `--fgf`: Serum FGF concentration in nM

## Reproducing Published Examples

To run the first replicate of cell population simulation with 0.003162 μM dose of trametinib for 72 hours, with 100 starting cells, using
the name 'in_silico_drs' and 16 CPUs, the following command can be entered:
mpirun -np 16 python cellpop_drs.py --sim_name in_silico_drs --cellpop 100 --exp_time 72 --drug trame_EC --dose 0.003162 --rep rep1
Upon completion of simulations, the results are saved to disk in a folder structure corresponding to drug name, replicate identifier and
drug dose respectively. For a single simulation with a specific replicate of a drug dose, outputs (temporal species trajectories)
from all cells in each generation are
saved in a python pickle object.

* To generate all drug dose response simulation data, run cellpop_drs.py for:

  * 4 drugs (alpelisib, neratinib, trametinib, palbociclib)
  * 10 dose levels (including control)
  * 10 replicates of each dose
* To generate cell population dynamics (number of alive cells over time) from dose response simulation outputs, run analysis_popdyn.py
  For this, results from all drug dose response simulations in step 5 need to be placed in the "output" folder in the main directory.
  Alternatively, outputs may be placed at a secondary locations and the path must be updated in line 68 of analysis_popdyn.py script.
  Outputs for the cell population dynamics will be saved in the "in_silico_drs_summary" folder under the output directory.
* To calculate GR score from the cell population dynamics, input files must be prepared for the gr-score calculation pipeline. For this, step 6
  must be completed first. After this, run analysis_grscore.py to generate the gr-score input file, which will be saved as 'drs_grcalc.tsv'
  in the "in_silico_drs_summary" folder.
* Take the input file generated at step 7 and run the gr-score calculation pipeline. For this:

  * 8a. Clone the gr-score git repository: https://github.com/datarail/gr_metrics
  * 8b. Install all dependencies including Python 2.0
  * 8c. Go to gr_metrics/SRC/python/scripts
  * 8d. Run python add_gr_column.py [input_path] > [output_path]
* Download all experimental dose response datasets (GR-scores) from here: https://www.synapse.org/#!Synapse:syn18456348/ and place them in
  'in_silico_drs_summary/mcf10a_drs_exp'
* Plots from figures 1,2 can be generated with jupyter notebooks included in the 'jupyter_notebooks' folder.

Running cell population simulation with a new single cell model:

By default, the cell population simulation workflow uses the SPARCED single cell model. It is capable of running simulations with a different single cell model given that the model has a compatible structure. A compatible model must satisfy the following requirements:

* The model must have a state matrix representing a single cell.
* The model must have a variable representing dynamic molecular signature of cell cycle markers, i.e., periodic activation and inactivation of cyclins.
* The model must be executable within a python module.

To replace the SPARCED model in cell population simulations with another single cell model:

1. Place all single cell simulation operations within a python function.
2. Write another python function to generate an input dict for the single cell model function, mirroring the input/output structure of the LoadSPARCED function.
3. Save both python functions as modules with the same name as the functions under bin/modules.
4. Write a json config file with key-specific values appropriate for the new model structure. Be sure to make "load_model" and "run_model" options consistent with the new module names. For more details on the stucture of the sim config, see sim_configs/README.md

The Novak-Tyson 1993 cell cycle model has been presented as an example for this procedure. The "load_model" and "run_model" modules have been provided as bin/modules/LoadTyson.py and bin/modules/RunTyson.py. The sim_config json file corresponding to this workflow is sim_config/default.json

## Contributors Guide

To enable broader portability of the LinResSims project, source code and dependencies are packaged into a distributable wheel using the `pyproject.toml` file and python's `build` command. Further, packages are installed at `/usr/local/lib/python3.10/site-packages/` via the `pip` package manager. In the event that you wish to contribute to update or change python packages, the distributable files (located at `LinResSims/dist/`) will need to be updated as well for the changes to take affect.

1. First, update the `pyproject.toml` file with any modifications to the package lists under the `dependencies` variable (line 18)
2. From the project root directory, execute the following command
   * `python -m build`
3. Lastly, install the new packages using pip:
   * `pipinstalldist/linressim-1.0-py3-none-any.whl --verbose --force`
