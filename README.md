<<<<<<< HEAD
# MD Setup for Carbonyl Probe Simulations
=======
# README: MD Setup for Carbonyl Probe Simulations
>>>>>>> d323832 (updates)

## Introduction

This repository provides a suite of scripts to automate the preparation of ligand topologies, conversion to GROMACS format, and the setup and analysis of Molecular Dynamics (MD) simulations in multiple solvents. The overall methodology is inspired by:

> Schneider, S. H.; Boxer, S. G. *Vibrational Stark Effects of Carbonyl Probes Applied to Reinterpret IR and Raman Data for Enzyme Inhibitors in Terms of Electric Fields at the Active Site.* **J. Phys. Chem. B** **2016**, *120* (36), 9672–9684.  
> DOI: [10.1021/acs.jpcb.6b08133](https://doi.org/10.1021/acs.jpcb.6b08133)

Additional elements are adapted from:
[https://github.com/yasminshamsudin/boxer-solvatochrom.git](https://github.com/yasminshamsudin/boxer-solvatochrom.git)

For solvent parameters, please refer to:

> van der Spoel, D.; Ghahremanpour, M. M.; Lemkul, J. *Small Molecule Thermochemistry: A Tool For Empirical Force Field Development*, **J. Phys. Chem. A** **2018**, *122*, 8982–8988.  
> Accessible via [Virtual Chemistry](https://virtualchemistry.org/ff.php#refs).

---

## Workflow Overview

The MD setup and analysis workflow is divided into several sequential steps:

1. **Topology Generation and Conversion (Step 1)**  
   **`step_1.sh`**  
   - Reads `smiles.txt` for each ligand and generates a 3D structure using Obabel.  
   - Assigns GAFF atom types via antechamber and parmchk2, and creates AMBER topologies using tleap.  
   - Converts AMBER files to GROMACS format with acpype.  
   - Cleans up and organizes outputs into a `step_1_outputs/` folder.

2. **MD Simulation Setup (Step 2)**  
   **`step_2.sh`**  
   - Uses the GROMACS files from Step 1 to build a simulation box and solvate the system in various solvents (H2O, D2O, and organic solvents).  
   - Generates both standard and zero-charge (0q) topologies for vibrational Stark analysis.  
   - Sets up triplicate MD simulation replicates (`run1/`, `run2/`, `run3/`) within each ligand’s solvent-specific MD folder, copying analysis scripts and creating an HPC submission script (`run_md.sh`).

3. **Job Submission (Step 3)**  
   **`step_3_modified.sh`**  
   - Iterates through each ligand’s MD run directories, ensures that the submission script (`run_md.sh`) is executable, changes to the appropriate directory, and submits the jobs via `sbatch`.  
   - This script ensures the working directory is correctly set so that all input files are found by GROMACS.

4. **Post-Simulation Analysis (Step 4)**  
   **`step_4.sh`**  
   - Loads the `scipy-stack` module to provide the necessary Python scientific libraries.  
   - Iterates through all run directories, executes `calc_fields.py` (which computes field data and generates `FIELDS.txt`), and then runs `time_series.py` with `FIELDS.txt` as input.  
   - Saves the output of `time_series.py` to `results.txt` in each run directory.

5. **Data Aggregation and Plotting**  
   Two Python scripts streamline the analysis of vibrational frequency data:  
   - **`aggregate_results.py`**  
     - Recursively parses each `results.txt` file (extracting key metrics from the "Eproj (average overall electric field)" block) from the MD runs, records one row per run, and then appends an average row per ligand–solvent pair.  
     - Saves the aggregated data to `aggregated_results.csv`.
   - **`plot_field_vs_frequency.py`**  
     - Reads `aggregated_results.csv` and frequency data from `freq_data.txt`.  
     - Merges the average electric field data with vibrational frequencies (provided per ligand and solvent).  
     - Creates a publication-quality scatter plot where the x-axis is the average solvent electric field (in MV/cm) and the y-axis is the C=O frequency (in cm⁻¹).  
     - For each ligand, a solid black regression line (extended across the full x-range) is fitted and printed to the terminal in the format:  
       ```
       v̅C═O<ligand> = <slope>|F⃗solv| + <intercept> (R2 = <r²>)
       ```  
     - The solvent markers are color-coded; a legend (showing only solvent names and colors) is placed outside to the right.  
     - The final plot is saved as `electric_field_vs_frequency.png`.

---

## Repository Structure

```
Auto-VSE/
├── ANALYSIS_SCRIPTS/
│   └── (custom analysis scripts, if any)
├── LIGANDS/
│   ├── ligandA/
│   │   ├── smiles.txt
│   │   ├── probe.ndx
│   │   └── (outputs from step_1.sh, step_2.sh, including *_md/ folders)
│   └── ligandB/
│       ├── smiles.txt
│       ├── probe.ndx
│       └── ...
├── MD_PARAMS/
│   ├── min.mdp
│   ├── nvt.mdp
│   ├── npt.mdp
│   └── md.mdp
├── SOLVENTS/
│   ├── H2O/
│   │   ├── spc216.gro
│   │   └── tip3p.itp
│   ├── D2O/
│   │   ├── spc216.gro
│   │   └── tip3p.itp (modified for D2O)
│   ├── dibutyl-ether/
│   │   ├── conf.gro
│   │   └── topol.top
│   └── (other solvents)/
│       ├── conf.gro
│       └── topol.top
├── freq_data.txt
├── step_1.sh
├── step_2.sh
├── step_3_modified.sh
├── step_4.sh
├── aggregate_results.py
└── plot_field_vs_frequency.py
```

**Key Files:**

<<<<<<< HEAD
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
=======
- **`step_1.sh`**: Converts SMILES to force-field parameterized topologies and prepares GROMACS input files.
- **`step_2.sh`**: Sets up triplicate MD simulations per ligand/solvent, prepares topologies (including zero-charge versions), and creates HPC submission scripts.
- **`step_3.sh`**: Iterates through simulation run directories, ensures executable permissions, and submits jobs via SLURM.
- **`step_4.sh`**: Loads the scipy-stack, executes `calc_fields.py` and `time_series.py` in each run directory, saving analysis output to `results.txt`.
- **`aggregate_results.py`**: Aggregates key vibrational frequency metrics from `results.txt` files into `aggregated_results.csv`.
- **`plot_field_vs_frequency.py`**: Merges electric field data with vibrational frequency data, fits regression lines per ligand, prints the regression equations, and creates a publication-quality plot (`electric_field_vs_frequency.png`).
>>>>>>> d323832 (updates)

---

## How to Run

<<<<<<< HEAD
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
=======
1. **Clone/Download** the repository into your working directory
  - git clone https://github.com/DomFico/Auto-VSE.git

2. **Prepare Ligands**:  
   - Create subdirectories under `LIGANDS/` for each ligand.  
   - Place `smiles.txt` (with one SMILES string on the first line) and `probe.ndx` in each ligand folder.

3. **Prepare Solvents**:  
   - Under `SOLVENTS/`, create a subdirectory for each solvent (e.g., H2O, D2O, trichloromethane, etc.).  
   - Place `spc216.gro`/`tip3p.itp` for water (or D2O) and `conf.gro`/`topol.top` for organic solvents in the corresponding folders.

4. **MD Parameters**:  
   - Ensure `MD_PARAMS/` contains the required `.mdp` files: `min.mdp`, `nvt.mdp`, `npt.mdp`, and `md.mdp`.

5. **Run the Workflow**:  
   - **Step 1**:  
>>>>>>> d323832 (updates)
     ```bash
     ./step_1.sh
     ```  
   - **Step 2**:  
     ```bash
     ./step_2.sh
     ```  
   - **Step 3**:  
     ```bash
     ./step_3.sh
     ```  
   - **Step 4**:  
     ```bash
     ./step_4.sh
     ```

6. **Data Aggregation and Plotting**:  
   - Aggregate the analysis results:  
     ```bash
     python aggregate_results.py
     ```  
   - Generate the publication-quality plot and print regression lines:  
     ```bash
     python plot_field_vs_frequency.py
     ```

---

## Required Software

- **Obabel** (for 3D geometry generation)  
- **AmberTools** (antechamber, parmchk2, tleap)  
- **Acpype** (for AMBER-to-GROMACS conversion)  
- **GROMACS**  
- **GNU Bash** and standard GNU utilities  
- **SLURM** (or equivalent HPC job scheduler)  
- **Python 3** with packages: `numpy`, `pandas`, `matplotlib`

---

## Final Notes

- **D2O Simulations**:  
  For D2O, the hydrogen mass is increased to 2.014 in the water model as an approximation.

- **HPC Job Submission**:  
  The SLURM directives in `run_md.sh` and `step_3.sh` may need adaptation to your cluster's specifications.

<<<<<<< HEAD
Have Fun :)
=======
- **Analysis Customization**:  
  Additional analysis scripts placed in `ANALYSIS_SCRIPTS/` are automatically copied into each run folder during Step 2.

For any issues or feedback, please open an issue or contact the maintainers.

Happy Simulating :)
```
>>>>>>> d323832 (updates)
