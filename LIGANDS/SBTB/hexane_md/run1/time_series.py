#!/usr/bin/env python3
"""
Script to compute autocorrelation-based error estimates.
It provides functions for computing the autocorrelation of a dataset,
combining autocorrelations from data halves, and estimating the error.
It also includes a random test for normally distributed data.

The expected data file (if it has 6 columns) is assumed to have the following columns:
  Column 0: Time (simulation time)
  Column 1: Eproj_C (electric field from the carbon atom)
  Column 2: Eproj_O (electric field from the oxygen atom)
  Column 3: Eproj (average overall electric field)
  Column 4: EprojDrop (electric field drop: difference between oxygen and carbon)
  Column 5: COlen (C–O bond length)
"""

import numpy as np
import random
import sys

def autocorr_1(x, npts):
    """
    Compute the autocorrelation function for array x with npts points.
    Note: Subtracting the mean introduces a slight anticorrelation.
    """
    m = np.mean(x)
    dx = x.copy()
    dx -= m
    result = np.zeros(npts)
    for i in range(npts):
        if i == 0:
            result[i] = np.mean(dx * dx)
        else:
            result[i] = np.mean(dx[i:] * dx[:-i])
        if (i + 1) % 10 == 0:
            # Print progress every 10 points (overwriting the same line)
            print("\rCorrelating point %i of %i : % .3e" % (i + 1, npts, result[i]), end='')
    print()  # Newline after progress output
    return result / result[0]

def autocorr_chop(x):
    """
    Compute a modified autocorrelation by combining the full dataset and its halves.
    """
    npts = min(10000, len(x) // 4)
    N = len(x)
    A0 = autocorr_1(x, npts)
    A1 = autocorr_1(x[:N // 2], npts)
    A2 = autocorr_1(x[N // 2:], npts)
    return 2 * A0 - 0.5 * (A1 + A2)

def chop(arr):
    """
    Returns a copy of the last 90% of the array.
    """
    return arr[-(9 * len(arr) // 10):].copy()

def est_err(data):
    """
    Estimate the error in the data using the autocorrelation function.
    Returns a tuple (Error Estimate, Gamma2 factor).
    """
    AC = autocorr_chop(data)
    Gam2 = 1.0
    for k, ack in enumerate(AC):
        if k == 0:
            continue
        Gam2 += 2 * (1 - float(k) / len(AC)) * ack
    if Gam2 < 1.0:  # Gamma2 must not be less than 1.
        Gam2 = 1.0
    ErrEst = Gam2 * np.std(data)**2 / len(data)
    return ErrEst, Gam2

def random_test():
    """
    Test the error estimation on data drawn from a normal distribution.
    """
    Samps = 100000
    # Using SystemRandom for better randomness
    r = random.SystemRandom()
    for npts in [10]:
        print("--- Testing for normal distribution, %i data points... ---" % npts)
        Gams = []
        for k in range(Samps):
            # Update progress on the same line
            print("\rSample %i" % k, end='')
            Data = np.array([random.gauss(0, 1) for i in range(npts)])
            _, Gam2 = est_err(Data)
            Gams.append(Gam2)
        Gams = np.array(Gams)
        print("\nAverage of Gam2 values:", np.mean(Gams))
        sortedGams = sorted(Gams)
        # Compute the interquartile range (IQR)
        IQR = sortedGams[3 * Samps // 4] - sortedGams[Samps // 4]
        FDC = 2 * IQR * Samps ** (-0.333)
        bins = int((max(Gams) - min(Gams)) / FDC)
        print("Calculated bins for histogram:", bins)
        # Uncomment below to plot histogram using matplotlib
        # import matplotlib.pyplot as plt
        # plt.hist(Gams, bins=bins)
        # plt.show()

def main():
    """
    Main function: loads data from file and computes error estimates for each column.
    """
    random_switch = 0
    if random_switch:
        random_test()
    else:
        if len(sys.argv) < 2:
            print("Usage: %s datafile" % sys.argv[0])
            sys.exit(1)
        # Load data from the provided file
        Data = np.loadtxt(sys.argv[1])
        
        # If the data file has 6 columns, print detailed column information.
        if Data.ndim == 2 and Data.shape[1] == 6:
            print("Data file columns detected: 6 columns")
            print("Column 0: Time (simulation time)")
            print("Column 1: Eproj_C (electric field at carbon)")
            print("Column 2: Eproj_O (electric field at oxygen)")
            print("Column 3: Eproj (average overall electric field)")
            print("Column 4: EprojDrop (difference between oxygen and carbon electric field)")
            print("Column 5: COlen (C-O bond length)")
            print("")
        else:
            print("Data file has %i columns. No detailed column descriptions provided." % (Data.shape[1] if Data.ndim == 2 else 1))
        
        # Process each column of the data file
        # If the file has multiple columns, compute error estimates for each.
        # Optionally, you could define a mapping of column index to description.
        descriptions = {
            0: "Time",
            1: "Eproj_C (electric field at carbon)",
            2: "Eproj_O (electric field at oxygen)",
            3: "Eproj (average overall electric field)",
            4: "EprojDrop (electric field drop: Eproj_O - Eproj_C)",
            5: "COlen (C-O bond length)"
        }
        
        if Data.ndim == 1:
            # If there is only one column, process it
            cols = [Data]
            col_indices = [0]
        else:
            cols = [Data[:, i] for i in range(Data.shape[1])]
            col_indices = range(Data.shape[1])
        
        for i in col_indices:
            # Print a header for this column
            if i in descriptions:
                print(f"Processing column {i}: {descriptions[i]}")
            else:
                print(f"Processing column {i}")
            for spac in [1]:
                print("Spacing =", spac)
                sliced_data = cols[i][::spac]
                print("Number of points =", len(sliced_data))
                ErrEst, Gam2 = est_err(sliced_data)
                print("Error enhancement factor (Gamma2) =", Gam2)
                print("Mean =", np.mean(sliced_data))
                print("Standard deviation =", np.std(sliced_data))
                print("Estimated error =", np.sqrt(ErrEst))
                print("")
                
if __name__ == "__main__":
    main()
