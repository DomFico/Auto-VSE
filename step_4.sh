#!/bin/bash
# -----------------------------------------------------------------------------
# step_4.sh
#
# This script loads the scipy-stack module, iterates through all MD run
# directories, executes calc_fields.py to compute field data (producing FIELDS.txt),
# then runs time_series.py with FIELDS.txt as input. The output (print statements)
# from time_series.py is saved as results.txt in each run directory.
#
# Usage: ./step_4.sh
#
# Note: Ensure that both calc_fields.py and time_series.py are present in each run
# directory and that sbatch or similar job schedulers are not interfering with file I/O.
# -----------------------------------------------------------------------------

# Load the scipy-stack module for the Python scientific environment
module load scipy-stack

# Get the directory where this script resides
script_dir=$(dirname "$(realpath "$0")")

# Define the LIGANDS directory based on the project structure
LIGANDS_DIR="$script_dir/LIGANDS"

# Check that the LIGANDS directory exists
if [[ ! -d "$LIGANDS_DIR" ]]; then
    echo "ERROR: LIGANDS directory not found at $LIGANDS_DIR"
    exit 1
fi

echo "Starting FIELDS calculation and time series analysis for all MD runs..."
echo "Ligands directory: $LIGANDS_DIR"

# Iterate over each ligand directory
for ligand_dir in "$LIGANDS_DIR"/*/; do
    ligand_name=$(basename "$ligand_dir")
    echo "Processing ligand: $ligand_name"

    # Iterate over each solvent MD folder (e.g., methanol_md/)
    for md_folder in "$ligand_dir"/*_md/; do
        if [[ ! -d "$md_folder" ]]; then
            echo "  No MD folder found in $ligand_name. Skipping."
            continue
        fi

        solvent_name=$(basename "$md_folder" _md)
        echo "  Processing solvent MD folder: ${solvent_name}_md"

        # Iterate over each run directory (e.g., run1, run2, run3)
        for run_dir in "$md_folder"/run*/; do
            if [[ ! -d "$run_dir" ]]; then
                echo "    No run directory found in ${solvent_name}_md. Skipping."
                continue
            fi
            run_name=$(basename "$run_dir")
            echo "    Processing run: $run_name for ligand: $ligand_name, solvent: $solvent_name"

            # Change into the run directory
            pushd "$run_dir" > /dev/null || { echo "      ERROR: Cannot change directory to $run_dir"; continue; }

            # Execute calc_fields.py if it exists
            if [[ -f "calc_fields.py" ]]; then
                echo "      Executing calc_fields.py..."
                python calc_fields.py
                if [[ $? -ne 0 ]]; then
                    echo "      ERROR: calc_fields.py execution failed in $(pwd)"
                    popd > /dev/null
                    continue
                fi
            else
                echo "      WARNING: calc_fields.py not found in $(pwd). Skipping this run."
                popd > /dev/null
                continue
            fi

            # Ensure that FIELDS.txt is produced
            if [[ ! -f "FIELDS.txt" ]]; then
                echo "      ERROR: FIELDS.txt not found after running calc_fields.py in $(pwd)"
                popd > /dev/null
                continue
            fi

            # Execute time_series.py using FIELDS.txt as input and redirect output to results.txt
            if [[ -f "time_series.py" ]]; then
                echo "      Executing time_series.py with FIELDS.txt; saving output to results.txt..."
                python time_series.py FIELDS.txt > results.txt 2>&1
                if [[ $? -eq 0 ]]; then
                    echo "      Time series analysis completed; results saved in results.txt"
                else
                    echo "      ERROR: time_series.py execution failed in $(pwd)"
                fi
            else
                echo "      WARNING: time_series.py not found in $(pwd). Skipping time series analysis."
            fi

            # Return to the previous directory
            popd > /dev/null
        done
    done
done

echo "FIELDS calculation and time series analysis completed for all MD runs."
