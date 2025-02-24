#!/usr/bin/env python3
"""
aggregate_results.py

This script navigates the LIGANDS directory structure, locates each run's results.txt file,
and extracts the following metrics from the block starting with:
    Processing column 3: Eproj (average overall electric field)
The metrics extracted are:
    Error enhancement factor (Gamma2) = <number>
    Mean = <number>
    Standard deviation = <number>
    Estimated error = <number>

For each ligand and solvent pair, the script writes one row per run
(with columns: Ligand, Solvent, Run, Gamma2, Mean, Std, EstError)
and then appends an additional row labeled "Average" that averages the three runs.
The aggregated data is saved as "aggregated_results.csv" in the current working directory.

Usage:
    python aggregate_results.py
"""

import os
import re
import csv
import statistics

# Define regex patterns for the markers
proc_pattern   = re.compile(r"Processing column\s+(\d+):\s+(.+)")
gamma_pattern  = re.compile(r"Error enhancement factor \(Gamma2\)\s*=\s*([\d\.eE\+\-]+)")
mean_pattern   = re.compile(r"Mean\s*=\s*([\d\.eE\+\-]+)")
std_pattern    = re.compile(r"Standard deviation\s*=\s*([\d\.eE\+\-]+)")
est_pattern    = re.compile(r"Estimated error\s*=\s*([\d\.eE\+\-]+)")

# We want to parse the block for column 3 with the following label:
desired_label = "Eproj (average overall electric field)"

# Dictionary to store results:
# results[ligand][solvent] = list of dicts, each dict holds metrics for a run with key "run"
results = {}

# Assume LIGANDS is in the current working directory
ligands_dir = os.path.join(os.getcwd(), "LIGANDS")
if not os.path.isdir(ligands_dir):
    print(f"ERROR: LIGANDS directory not found at {ligands_dir}")
    exit(1)

print(f"Scanning results under: {ligands_dir}")

# Traverse each ligand directory
for ligand in os.listdir(ligands_dir):
    ligand_path = os.path.join(ligands_dir, ligand)
    if not os.path.isdir(ligand_path):
        continue
    results.setdefault(ligand, {})
    # Solvent folders end with _md
    for folder in os.listdir(ligand_path):
        if not folder.endswith("_md"):
            continue
        solvent = folder[:-3]  # remove the trailing "_md"
        solvent_path = os.path.join(ligand_path, folder)
        if not os.path.isdir(solvent_path):
            continue
        results[ligand].setdefault(solvent, [])
        # Process each run folder (e.g., run1, run2, run3)
        for run in os.listdir(solvent_path):
            run_path = os.path.join(solvent_path, run)
            if not os.path.isdir(run_path):
                continue
            results_file = os.path.join(run_path, "results.txt")
            if not os.path.isfile(results_file):
                print(f"WARNING: results.txt not found in {run_path}. Skipping.")
                continue

            # Read and parse results.txt
            with open(results_file, "r") as f:
                lines = f.readlines()

            in_block = False
            metrics = {}
            for line in lines:
                line_strip = line.strip()
                proc_match = proc_pattern.match(line_strip)
                if proc_match:
                    col_num, label = proc_match.groups()
                    # If we are already in the desired block and encounter a new block, stop parsing
                    if in_block:
                        break
                    if desired_label in label:
                        in_block = True
                        continue
                if in_block:
                    if "gamma2" not in metrics:
                        gamma_match = gamma_pattern.search(line_strip)
                        if gamma_match:
                            metrics["gamma2"] = float(gamma_match.group(1))
                    if "mean" not in metrics:
                        mean_match = mean_pattern.search(line_strip)
                        if mean_match:
                            metrics["mean"] = float(mean_match.group(1))
                    if "std" not in metrics:
                        std_match = std_pattern.search(line_strip)
                        if std_match:
                            metrics["std"] = float(std_match.group(1))
                    if "est_error" not in metrics:
                        est_match = est_pattern.search(line_strip)
                        if est_match:
                            metrics["est_error"] = float(est_match.group(1))
            if in_block and len(metrics) == 4:
                print(f"Parsed {results_file} for ligand '{ligand}', solvent '{solvent}', run '{run}'.")
                # Store the run name along with metrics
                metrics["run"] = run
                results[ligand][solvent].append(metrics)
            else:
                print(f"WARNING: Desired metrics not found in {results_file}. Skipping this run.")

# Build rows for CSV output.
# Each row: Ligand, Solvent, Run, Gamma2, Mean, Std, EstError
rows = []
for ligand, solvent_dict in results.items():
    for solvent, runs in solvent_dict.items():
        # Add one row per run
        for run_data in runs:
            row = {
                "Ligand": ligand,
                "Solvent": solvent,
                "Run": run_data["run"],
                "Gamma2": run_data["gamma2"],
                "Mean": run_data["mean"],
                "Std": run_data["std"],
                "EstError": run_data["est_error"],
            }
            rows.append(row)
        # Add a row for the average across runs
        if runs:
            gamma2_avg = statistics.mean([r["gamma2"] for r in runs])
            mean_avg   = statistics.mean([r["mean"] for r in runs])
            std_avg    = statistics.mean([r["std"] for r in runs])
            est_avg    = statistics.mean([r["est_error"] for r in runs])
            avg_row = {
                "Ligand": ligand,
                "Solvent": solvent,
                "Run": "Average",
                "Gamma2": gamma2_avg,
                "Mean": mean_avg,
                "Std": std_avg,
                "EstError": est_avg,
            }
            rows.append(avg_row)

# Write aggregated data to CSV file
csv_filename = "aggregated_results.csv"
csv_fields = ["Ligand", "Solvent", "Run", "Gamma2", "Mean", "Std", "EstError"]

with open(csv_filename, "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_fields)
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

print(f"Aggregation complete. Results saved to {csv_filename}.")
