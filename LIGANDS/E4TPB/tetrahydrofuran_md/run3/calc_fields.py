#!/usr/bin/env python3
"""
Script to calculate electric fields for a solute in solvent with a carbonyl vibrational probe.
Reads coordinates and forces from .xvg files and writes the calculated fields to FIELDS.txt.
"""

import sys
import numpy as np

def isnumber(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

# Use context managers for safe file handling
try:
    with open('co_coords.xvg', 'r') as f_x, \
         open('co_forces.xvg', 'r') as f_f, \
         open('co_forces_0q.xvg', 'r') as f_f0q, \
         open('FIELDS.txt', 'w') as fields:
        
        for line_x, line_f, line_f0q in zip(f_x, f_f, f_f0q):
            tokens_x = line_x.split()
            # Skip non-data lines (e.g., comments) where the first token is not a number
            if tokens_x and isnumber(tokens_x[0]):
                # Get the time values from the coordinate and discharged forces files
                time_x = float(tokens_x[0])
                time_f0q = float(line_f0q.split()[0])
                if time_x != time_f0q:
                    print("Time indices do not match. Exiting.")
                    sys.exit(1)

                # Parse C and O coordinates; tokens 1-3 for C and 4-6 for O
                try:
                    x_C = np.array([float(val) for val in tokens_x[1:4]])
                    x_O = np.array([float(val) for val in tokens_x[4:7]])
                except Exception as e:
                    print("Error parsing coordinates:", e)
                    sys.exit(1)
                
                # Compute the CO bond vector, its length, and the unit vector
                COvec = x_O - x_C
                COlen = np.linalg.norm(COvec)
                if COlen == 0:
                    print("Zero bond length encountered. Exiting.")
                    sys.exit(1)
                COunitvec = COvec / COlen

                # Parse forces: tokens 1-3 for C and 4-6 for O in both force files
                try:
                    f_tokens = line_f.split()
                    f0q_tokens = line_f0q.split()
                    f_C = np.array([float(val) for val in f_tokens[1:4]])
                    f_O = np.array([float(val) for val in f_tokens[4:7]])
                    f0q_C = np.array([float(val) for val in f0q_tokens[1:4]])
                    f0q_O = np.array([float(val) for val in f0q_tokens[4:7]])
                except Exception as e:
                    print("Error parsing forces:", e)
                    sys.exit(1)

                # Calculate the net electrostatic forces by subtracting discharged forces
                fe_C = f_C - f0q_C
                fe_O = f_O - f0q_O

                # Project the forces onto the CO bond direction
                feproj_C = np.dot(fe_C, COunitvec)
                feproj_O = np.dot(fe_O, COunitvec)

                # Convert the force projections to electric fields using provided constants
                Eproj_C = (feproj_C / 0.589301) * 0.1036427
                Eproj_O = (feproj_O / -0.539101) * 0.1036427
                Eproj = (Eproj_C + Eproj_O) / 2
                EprojDrop = (Eproj_O - Eproj_C)

                # Prepare and write the output line to FIELDS.txt
                fieldInfo = f"{time_x}\t{Eproj_C}\t{Eproj_O}\t{Eproj}\t{EprojDrop}\n"
                fields.write(fieldInfo)
        
        print("Field calculations completed successfully.")
except IOError as e:
    print("File I/O error:", e)
    sys.exit(1)
