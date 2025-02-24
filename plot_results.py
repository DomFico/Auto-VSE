#!/usr/bin/env python3
"""
plot_field_vs_frequency.py

This script plots vibrational frequencies versus average electric field values.
It reads the aggregated_results.csv (with average electric field values per ligand–solvent pair)
and parses freq_data.txt (with vibrational frequencies for each ligand–solvent pair).
Only ligand–solvent pairs present in both files are plotted.

For each ligand:
    - Data points are shown as scatter markers colored by solvent.
    - A solid regression line (in black) is computed and extended across the full x-range.
    - The regression line is printed to the terminal in the format:
          v̅C═O<ligand> = <slope>|F⃗solv| + <intercept> (R2 = <r²>)
    
The x-axis is labeled:
    "Average Solvent Electric Field / (MV/cm)"
and the y-axis is labeled:
    "C=O Frequency / cm⁻¹"

A legend (without a title) showing the solvent marker colors is placed outside (to the right) of the graph.
The final plot is saved as a PNG file (electric_field_vs_frequency.png).

Usage:
    python plot_field_vs_frequency.py
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Function to parse frequency data from freq_data.txt ---
def parse_freq_data(filepath):
    """
    Parses a frequency file with the following structure:
    
    LIGAND
        Solvent: value±error
        ...
    
    Returns a dictionary: { ligand: { solvent: (frequency, error) } }
    """
    freq_data = {}
    current_ligand = None
    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            # Lines without leading whitespace denote a ligand name.
            if not line.startswith(' '):
                current_ligand = stripped
                freq_data[current_ligand] = {}
            else:
                # Expect format like: "D2O: 1703.8±0.3"
                m = re.match(r"(\S+):\s+(\S+)", stripped)
                if m:
                    solvent, value_str = m.groups()
                    if "N/A" in value_str:
                        continue  # Skip missing data
                    parts = value_str.split("±")
                    if len(parts) == 2:
                        try:
                            freq = float(parts[0])
                            err  = float(parts[1])
                            freq_data[current_ligand][solvent] = (freq, err)
                        except ValueError:
                            continue
    return freq_data

# --- Load aggregated electric field data ---
csv_filename = "aggregated_results.csv"
if not os.path.isfile(csv_filename):
    raise FileNotFoundError(f"{csv_filename} not found in the current directory.")

# Read the CSV and filter for rows where Run == "Average"
df_field = pd.read_csv(csv_filename)
df_field = df_field[df_field["Run"] == "Average"]

# --- Parse frequency data ---
freq_file = "freq_data.txt"
if not os.path.isfile(freq_file):
    raise FileNotFoundError(f"{freq_file} not found in the current directory.")
freq_data = parse_freq_data(freq_file)

# --- Merge datasets ---
# Build a dictionary: merged_data[ligand] = list of tuples (solvent, avg_field, frequency)
merged_data = {}
for _, row in df_field.iterrows():
    ligand = row["Ligand"]
    solvent = row["Solvent"]
    avg_field = row["Mean"]  # Average electric field value
    # Merge only if frequency data exists for this ligand and solvent
    if ligand in freq_data and solvent in freq_data[ligand]:
        frequency = freq_data[ligand][solvent][0]
        merged_data.setdefault(ligand, []).append((solvent, avg_field, frequency))

# --- Define color mappings ---
# Fixed solvent order for consistency
solvent_order = ["D2O", "trichloromethane", "dichloromethane",
                 "dimethyl-sulfoxide", "tetrahydrofuran", "dibutyl-ether", "hexane"]
cmap_solvent = plt.get_cmap("tab10")
solvent_colors = {sol: cmap_solvent(i % 10) for i, sol in enumerate(solvent_order)}

# --- Determine global x-range from all data points ---
all_x = []
for ligand, data in merged_data.items():
    for (sol, x_val, freq) in data:
        all_x.append(x_val)
if not all_x:
    raise ValueError("No data points found for merging.")
global_xmin, global_xmax = min(all_x), max(all_x)
x_range_plot = np.linspace(global_xmin, global_xmax, 200)

# --- Set up publication-quality parameters ---
plt.rcParams.update({
    "font.size": 12,
    "figure.figsize": (12, 6),  # Wider figure
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 10,
    "lines.linewidth": 2,
    "axes.linewidth": 1.5,
})
plt.figure()

# --- Plotting ---
# For each ligand, plot its data points and a regression line (solid black)
for ligand, data in merged_data.items():
    # Sort data by average electric field
    data.sort(key=lambda x: x[1])
    solvents_list, x_vals, y_vals = zip(*data)
    x_arr = np.array(x_vals)
    y_arr = np.array(y_vals)
    
    # Perform linear regression if at least 2 points exist
    if len(x_arr) >= 2:
        a, b = np.polyfit(x_arr, y_arr, 1)
        y_fit = a * x_range_plot + b
        ss_res = np.sum((y_arr - (a * x_arr + b)) ** 2)
        ss_tot = np.sum((y_arr - np.mean(y_arr)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot != 0 else 0
    else:
        a, b, r2 = 0, 0, 0

    # Print the regression line for this ligand to the terminal
    print(f"v̅C═O{ligand} = {a:.2f}|F⃗solv| + {b:.1f} (R2 = {r2:.2f})")
    
    # Plot each data point (scatter) colored by its solvent with black edge
    for sol, x_pt, y_pt in data:
        color_marker = solvent_colors.get(sol, "black")
        plt.scatter(x_pt, y_pt, color=color_marker, s=80, edgecolor="k", zorder=3)
    
    # Plot the regression line as a solid black line spanning the full x-range
    if len(x_arr) >= 2:
        plt.plot(x_range_plot, y_fit, linestyle="-", color="black", linewidth=1.5)

# --- Customize axes and labels ---
plt.xlabel("Average Solvent Electric Field / (MV/cm)")
plt.ylabel("C=O Frequency / cm⁻¹")
plt.xlim(global_xmin - 0.05*(global_xmax - global_xmin), global_xmax + 0.05*(global_xmax - global_xmin))
# Remove title

# --- Create legend for solvents (without a title) ---
legend_handles = []
used_solvents = set()
for ligand, data in merged_data.items():
    for sol, _, _ in data:
        used_solvents.add(sol)
used_solvents = [sol for sol in solvent_order if sol in used_solvents]
for sol in used_solvents:
    handle = plt.Line2D([], [], marker='o', linestyle='', color=solvent_colors[sol],
                          markersize=8, label=sol)
    legend_handles.append(handle)
plt.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)

# Remove grid
# plt.grid(False) is not needed as grid is off by default in this configuration

plt.tight_layout()

# Save the plot as a publication-quality PNG file
output_png = "electric_field_vs_frequency.png"
plt.savefig(output_png, dpi=300, bbox_inches="tight")
print(f"Plot saved as {output_png}")
plt.show()
