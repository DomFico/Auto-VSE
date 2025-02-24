#!/bin/bash
# -----------------------------------------------------------------------------
# step_3_modified.sh
#
# This script iterates through each ligand directory under the LIGANDS folder,
# then for each solvent MD folder (ending with _md) and for each run directory,
# it changes into that run folder, makes sure run_md.sh is executable, and then
# submits the MD job using sbatch.
#
# Changing into the run directory sets the working directory for the job correctly,
# ensuring that all necessary input files (e.g., mdp, gro, top) are accessible.
#
# Usage: ./step_3_modified.sh
#
# Note: Ensure that sbatch is available in your environment.
# -----------------------------------------------------------------------------

# Get the directory where this script resides
script_dir=$(dirname "$(realpath "$0")")

# Define the LIGANDS directory
LIGANDS_DIR="$script_dir/LIGANDS"

# Check that LIGANDS directory exists
if [[ ! -d "$LIGANDS_DIR" ]]; then
    echo "ERROR: LIGANDS directory not found at $LIGANDS_DIR"
    exit 1
fi

echo "Starting job submission with correct working directories..."
echo "Ligands directory: $LIGANDS_DIR"

# Iterate over each ligand directory
for ligand_dir in "$LIGANDS_DIR"/*/; do
    ligand_name=$(basename "$ligand_dir")
    echo "--------------------------------------------"
    echo "Processing ligand: $ligand_name"

    # Iterate over each solvent MD folder (e.g., methanol_md/)
    for md_folder in "$ligand_dir"/*_md/; do
        if [[ ! -d "$md_folder" ]]; then
            echo "  No MD folder found in $ligand_dir. Skipping."
            continue
        fi

        # Extract solvent name by removing the trailing '_md'
        solvent_name=$(basename "$md_folder" _md)
        echo "  Processing solvent MD folder: ${solvent_name}_md"

        # Iterate over each run directory (run1, run2, run3)
        for run_dir in "$md_folder"/run*/; do
            if [[ ! -d "$run_dir" ]]; then
                echo "    No run directory found in $md_folder. Skipping."
                continue
            fi
            run_name=$(basename "$run_dir")
            echo "    Submitting job for: Ligand: $ligand_name | Solvent: ${solvent_name} | Run: $run_name"

            # Change into the run directory to set the working directory correctly
            pushd "$run_dir" > /dev/null || { echo "      ERROR: Cannot change directory to $run_dir"; continue; }

            # Set executable permissions for run_md.sh
            if [[ -f "run_md.sh" ]]; then
                chmod +x run_md.sh
            fi

            # Check if the submission script exists and is executable
            if [[ -x "run_md.sh" ]]; then
                echo "      Submitting job from $(pwd)"
                sbatch_output=$(sbatch run_md.sh 2>&1)
                if [[ $? -eq 0 ]]; then
                    echo "      Job submitted successfully: $sbatch_output"
                else
                    echo "      ERROR: Job submission failed in $(pwd)"
                    echo "      sbatch output: $sbatch_output"
                fi
            else
                echo "      WARNING: run_md.sh not found or not executable in $(pwd). Skipping this run."
            fi

            # Return to the previous directory
            popd > /dev/null
        done
    done
done

echo "Job submission process completed for all MD runs."
