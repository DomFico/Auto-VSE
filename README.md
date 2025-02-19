Here’s the corrected markdown that will render properly on GitHub:

```markdown
# README: MD Setup for Carbonyl Probe Simulations

## Introduction

This repository provides a set of scripts to automate the process of preparing ligand topologies, converting them to GROMACS format, and setting up Molecular Dynamics (MD) simulations in multiple solvents. These scripts take inspiration from the methodology discussed in the following work:

> Schneider, S. H.; Boxer, S. G. *Vibrational Stark Effects of Carbonyl Probes Applied to Reinterpret IR and Raman Data for Enzyme Inhibitors in Terms of Electric Fields at the Active Site.* **J. Phys. Chem. B** **2016**, *120* (36), 9672–9684.  
> DOI: [10.1021/acs.jpcb.6b08133](https://doi.org/10.1021/acs.jpcb.6b08133)

Additionally, portions of the setup process and structural files are adapted from:  
[https://github.com/yasminshamsudin/boxer-solvatochrom.git](https://github.com/yasminshamsudin/boxer-solvatochrom.git)

For the solvent parameter files, we rely on data described in:

> David van der Spoel, Mohammad M. Ghahremanpour, and Justin Lemkul,  
> *Small Molecule Thermochemistry: A Tool For Empirical Force Field Development*,  
> **J. Phys. Chem. A** **2018**, *122*, 8982–8988.  

These parameters can be accessed at [Virtual Chemistry](https://virtualchemistry.org/ff.php#refs). If you wish to use any of those solvent force fields, simply download the corresponding files (`topol.top`, `conf.gro`) and place them into your `SOLVENTS/` directory, creating a subfolder for each solvent.

---

## Overview of the Workflow

The MD workflow is split into two main scripts:

1. **`step_1.sh`**  
   - Converts SMILES strings into force-field parameterized topologies (`.top`, `.gro`) suitable for GROMACS.  
   - Uses `obabel`, `antechamber`, and `parmchk2` to assign GAFF atom types.  
   - Calls `tleap` (AmberTools) to generate AMBER topology files, then `acpype` to convert them into GROMACS files.  
   - Cleans intermediate files and organizes outputs into a `step_1_outputs` directory.

2. **`step_2.sh`**  
   - Sets up triplicate MD simulations for each ligand in multiple solvents (e.g., H2O, D2O, organic solvents).  
   - Places each simulation into dedicated folders (e.g., `solvent_md/run1`, `run2`, `run3`) within each ligand directory.  
   - Modifies topology files to handle zero-charge (0q) topologies, especially for calculating vibrational Stark effects.  
   - Copies analysis scripts to each run directory and creates a submission script (`run_md.sh`) for HPC job submission.

---

## Repository Structure

Below is an outline of the directory hierarchy expected by the scripts:

```
carbonyl_probe/
├── ANALYSIS_SCRIPTS/
│   └── (... your custom analysis scripts ...)
├── LIGANDS/
│   ├── ligandA/
│   │   ├── smiles.txt
│   │   ├── probe.ndx
│   │   └── (will contain outputs from step_1.sh and step_2.sh)
│   └── ligandB/
│       ├── smiles.txt
│       ├── probe.ndx
│       └── ...
├── MD_PARAMS/
│   ├── min.mdp
│   ├── md.mdp
│   ├── npt.mdp
│   └── nvt.mdp
├── SOLVENTS/
│   ├── H2O/
│   │   ├── spc216.gro
│   │   └── tip3p.itp
│   ├── D2O/
│   │   ├── spc216.gro (or custom D2O structure)
│   │   └── tip3p.itp (modified for D2O)
│   ├── dibutyl-ether/
│   │   ├── conf.gro
│   │   └── topol.top
│   └── (other solvents)/
│       ├── conf.gro
│       └── topol.top
├── step_1.sh
└── step_2.sh
```

**Important files and directories**:

- **`LIGANDS/`**  
  Each ligand must be placed in its own subdirectory. The minimal requirements inside each ligand folder are:
  - `smiles.txt` containing the SMILES string of the ligand (line 1).  
  - `probe.ndx` containing an index group file for the carbonyl probe atoms (e.g., for post-MD analysis).  
  - Any outputs from `step_1.sh` and `step_2.sh` will be placed here (e.g., `*_GMX.gro`, `*_GMX.top`, `*_md/` folders, etc.).

- **`MD_PARAMS/`**  
  Contains the GROMACS `.mdp` files for each MD phase: energy minimization (`min.mdp`), NVT equilibration (`nvt.mdp`), NPT equilibration (`npt.mdp`), and production MD (`md.mdp`).

- **`SOLVENTS/`**  
  Each solvent is placed in its own folder with at least:
  - `conf.gro` and `topol.top` for non-water solvents.  
  - `spc216.gro` and `tip3p.itp` for water (`H2O`), or similar for `D2O`.  
  The scripts use these to solvate the ligand and build the final topologies.  
  If you need additional solvent parameter files, you can obtain them from [Virtual Chemistry](https://virtualchemistry.org/ff.php#refs) and place them here.

- **`ANALYSIS_SCRIPTS/`** (Optional but recommended)  
  Place any custom analysis scripts here. They are automatically copied into every `run*` folder created by `step_2.sh`.

- **`step_1.sh`**  
  Processes each ligand:  
  1. Reads `smiles.txt` → generates a `.mol2` file (via `obabel`).  
  2. Assigns GAFF atom types (via `antechamber` + `parmchk2`).  
  3. Creates AMBER topologies (`tleap`) → converts to GROMACS (`acpype`).  
  4. Cleans up intermediate files into `step_1_outputs/`.

- **`step_2.sh`**  
  Uses the final `.gro` and `.top` from Step 1 to:  
  1. Create a simulation box around each ligand.  
  2. Insert solvent molecules or water boxes.  
  3. Prepare zero-charge (0q) topologies for vibrational Stark analysis.  
  4. Split the simulations into three replicates (`run1`, `run2`, `run3`).  
  5. Copy the analysis scripts from `ANALYSIS_SCRIPTS/`.  
  6. Generate a submission script (`run_md.sh`) for HPC job execution.

---

## Special Note on D2O

For D2O simulations included here, we simply increase the hydrogen mass to 2.014 within the water model. This approach is a rough approximation and not strictly correct, but it may provide a first pass at simulating heavy water in these simulations.

---

## How to Run

1. **Clone or Download** this repository into a working directory (e.g., `carbonyl_probe/`).  
2. **Place Your Ligands**:  
   - Create subfolders in `LIGANDS/`, one for each ligand.  
   - Put `smiles.txt` (with exactly one SMILES on the first line) in each ligand folder.  
   - Include `probe.ndx` to specify atomic indices for the carbonyl probe.  
3. **Check Solvents**:  
   - In `SOLVENTS/`, create subfolders for each solvent.  
   - Ensure that each subfolder has either `spc216.gro` + `tip3p.itp` (for H2O/D2O) or `conf.gro` + `topol.top` (for non-water).  
   - You may drag and drop additional solvent parameter files from [Virtual Chemistry](https://virtualchemistry.org/ff.php#refs) if needed.  
4. **Check MD Parameters**:  
   - Ensure that `MD_PARAMS/` contains `.mdp` files: `min.mdp`, `nvt.mdp`, `npt.mdp`, and `md.mdp`.  
5. **Run Step 1**:  
   - From the `carbonyl_probe/` directory, execute:  
     ```bash
     ./step_1.sh
     ```  
   - This will process all ligands in `LIGANDS/`. You will end up with `.gro` and `.top` files for each ligand, plus a `step_1_outputs/` folder with logs and intermediate files.  
6. **Run Step 2**:  
   - From the same `carbonyl_probe/` directory, run:  
     ```bash
     ./step_2.sh
     ```  
   - This will set up triplicate MD runs for each ligand in each solvent, creating a `solvent_name_md/` folder under each ligand. Each replicate run is inside `run1/`, `run2/`, `run3/`.  

---

## Required Software

- **Obabel**: 3D geometry generation from SMILES.  
- **AmberTools** (including `antechamber`, `parmchk2`, and `tleap`).  
- **Acpype**: Convert AMBER topologies to GROMACS format.  
- **GROMACS**: For box setup, solvation, and MD simulations.  
- **GNU Bash** (and standard GNU utilities: `grep`, `awk`, `sed`, etc.).

---

## Final Notes

- **Job Submission**:  
  Each `run_md.sh` script is configured for SLURM with a single GPU and 4 CPU cores. Adapt the `#SBATCH` headers and `module load` commands as needed for your HPC environment.  
- **Analysis**:  
  After each production MD, the script `run_md.sh` (in each replicate folder) reruns the trajectory on the zero-charge (0q) topology to collect vibrational Stark effect data (forces on the carbonyl probe). Any additional scripts placed in `ANALYSIS_SCRIPTS/` will be copied into every run folder automatically.  
- **Modifications**:  
  You may adapt the `.mdp` files in `MD_PARAMS/` or the HPC directives in `run_md.sh` as your system and queueing requirements dictate.  

Feel free to open an issue or contact the authors if you encounter any problems or if you have feedback on making these scripts more robust.

Have Fun :)
```
