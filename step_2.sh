#!/bin/bash
# -----------------------------------------------------------------------------
# step_2.sh
#
# This script prepares triplicate MD simulations for each ligand in multiple
# solvents. For each ligand/solvent combination, it:
#   1) Creates an MD folder named <solvent_name>_md/ containing subfolders
#      run1/, run2/, and run3/.
#   2) Copies over solvent configuration (e.g., conf.gro/topol.top or spc216.gro/tip3p.itp).
#   3) Generates zero-charge (0q) topologies by setting atomic charges to 0.0
#      in the new .itp files or tip3p_0q.itp (for water).
#   4) Solvates the system, modifies topologies (splitting into standard and
#      0q versions).
#   5) Copies any analysis scripts to each run folder.
#   6) Creates a submission script (run_md.sh) for each run with suitable
#      MD steps (min, nvt, npt, md, plus zero-charge rerun).
#
# Directory structure (assumed):
#
# carbonyl_probe/
# ├── LIGANDS
# │      ├── <ligand1>/
# │      │      ├── <ligand1>_GMX.gro
# │      │      ├── <ligand1>_GMX.top
# │      │      ├── probe.ndx
# │      │      └── ... (other files)
# │      └── <ligand2>/ ...
# ├── MD_PARAMS
# │      ├── min.mdp
# │      ├── npt.mdp
# │      ├── nvt.mdp
# │      └── md.mdp
# ├── SOLVENTS
# │      ├── H2O/
# │      │      ├── spc216.gro
# │      │      └── tip3p.itp
# │      ├── D2O/
# │      │      ├── spc216.gro
# │      │      └── tip3p.itp
# │      └── dibutyl-ether/ (other organic solvent examples)
# │             ├── conf.gro
# │             └── topol.top
# ├── ANALYSIS_SCRIPTS
# │      └── (any analysis scripts you want copied to each run folder)
# ├── step_1.sh
# └── step_2.sh
# -----------------------------------------------------------------------------

# Capture the directory that this script resides in
script_dir=$(dirname "$(realpath "$0")")

# Define key paths
LIGANDS_DIR="$script_dir/LIGANDS"
SOLVENTS_DIR="$script_dir/SOLVENTS"
MD_PARAMS_DIR="$script_dir/MD_PARAMS"
ANALYSIS_SCRIPTS_DIR="$script_dir/ANALYSIS_SCRIPTS"

# Check if the ANALYSIS_SCRIPTS directory exists; this is optional
if [[ ! -d "$ANALYSIS_SCRIPTS_DIR" ]]; then
    echo "WARNING: ANALYSIS_SCRIPTS directory not found at $ANALYSIS_SCRIPTS_DIR"
    echo "Will continue without copying analysis scripts"
else
    echo "Found ANALYSIS_SCRIPTS directory at $ANALYSIS_SCRIPTS_DIR"
    file_count=$(find "$ANALYSIS_SCRIPTS_DIR" -type f | wc -l)
    echo "Found $file_count files in ANALYSIS_SCRIPTS directory"
fi

# Iterate over each ligand subdirectory in LIGANDS
for ligand_dir in "$LIGANDS_DIR"/*/; do
    ligand_name=$(basename "$ligand_dir")
    echo "=============================================="
    echo "Processing ligand: $ligand_name"

    # Change to the ligand directory; skip if it fails
    cd "$ligand_dir" || { echo "ERROR: Cannot cd into $ligand_dir"; continue; }
    
    # Remove any previous *_md directories to allow a fresh setup
    rm -rf *_md

    # Check that the ligand's GROMACS files exist
    if [[ ! -f "${ligand_name}_GMX.gro" ]]; then
        echo "ERROR: ${ligand_name}_GMX.gro not found. Skipping ligand."
        cd "$script_dir" || exit
        continue
    fi

    # Create a clean backup of the original topology file
    cp "${ligand_name}_GMX.top" "${ligand_name}_GMX_original.top"
    echo "Created clean backup of original topology file."

    # Capture the absolute path of the ligand directory for reference
    ligand_abs=$(pwd)
    
    # Iterate over each solvent directory in SOLVENTS
    for solvent_dir in "$SOLVENTS_DIR"/*/; do
        solvent_name=$(basename "$solvent_dir")
        echo "  Processing solvent: $solvent_name for ligand: $ligand_name"
        
        # Create the parent MD folder for this solvent, e.g., methanol_md/
        parent_md_folder="${ligand_abs}/${solvent_name}_md"
        mkdir -p "$parent_md_folder"
        echo "    Created parent MD folder: $parent_md_folder"
        
        # Create triplicate runs: run1, run2, run3
        for run_num in 1 2 3; do
            run_folder="${parent_md_folder}/run${run_num}"
            mkdir -p "$run_folder"
            echo "      Created run folder: $run_folder"
            
            # Copy the appropriate solvent configuration files into the run folder
            if [[ "$solvent_name" == "H2O" || "$solvent_name" == "D2O" ]]; then
                cp "$solvent_dir/spc216.gro" "$run_folder/"
                cp "$solvent_dir/tip3p.itp" "$run_folder/"
                echo "      Copied H2O/D2O files to run${run_num}."
            else
                cp "$solvent_dir/conf.gro" "$run_folder/"
                cp "$solvent_dir/topol.top" "$run_folder/"
                echo "      Copied solvent conf.gro and topol.top to run${run_num}."
            fi

            # Enter the run folder
            cd "$run_folder" || { echo "      ERROR: Cannot cd into $run_folder"; continue; }
            
            # For water solvents, create tip3p_0q.itp by copying tip3p.itp and setting charges to 0.0
            if [[ "$solvent_name" == "H2O" || "$solvent_name" == "D2O" ]]; then
                cp "$solvent_dir/tip3p.itp" tip3p_0q.itp
                echo "      Created tip3p_0q.itp with zero charges for water."
                # Modify the [ atoms ] section in tip3p_0q.itp to set the charge field (column 7) to 0.0
                awk 'BEGIN { in_atoms=0 }
                     /^\[ *atoms *\]/ { in_atoms=1; print; next }
                     in_atoms && /^\[/ { in_atoms=0 }
                     in_atoms && NF > 0 {
                         $7="0.0";  # Set the charge to 0.0
                         print;
                         next
                     }
                     { print }' tip3p_0q.itp > tmp_tip3p_0q.itp && mv tmp_tip3p_0q.itp tip3p_0q.itp
            fi
            
            # Step 2: Create a cubic simulation box around the ligand with 2.0 nm padding
            editconf_cmd="gmx editconf -f ../../${ligand_name}_GMX.gro -o ${ligand_name}.gro -bt cubic -d 2.0 -c"
            echo "      Running: $editconf_cmd"
            eval $editconf_cmd

            # Always start with a clean copy of the original topology for solvation
            cp "../../${ligand_name}_GMX_original.top" "./temp_topol.top"

            # Step 3: Solvate the box
            if [[ "$solvent_name" == "H2O" || "$solvent_name" == "D2O" ]]; then
                # For water, use gmx solvate with spc216.gro and track the updated temp_topol.top
                solvate_cmd="gmx solvate -cp ${ligand_name}.gro -cs spc216.gro -o ${ligand_name}_${solvent_name}.gro -p temp_topol.top"
                echo "      Running: $solvate_cmd"
                $solvate_cmd 2>&1 | tee solvate.log
                
                # Parse the number of water molecules from the updated topology
                solvent_mol=$(grep -A 10 '^\[ molecules \]' temp_topol.top | grep 'SOL' | awk '{print $2}')
                if [[ -z "$solvent_mol" ]]; then
                    echo "      WARNING: Could not parse SOL molecule count. Defaulting to 1000."
                    solvent_mol=1000
                fi
            else
                # For non-water solvents, use gmx insert-molecules with conf.gro as the solvent
                solvate_cmd="gmx insert-molecules -f ${ligand_name}.gro -ci conf.gro -nmol 1000 -o ${ligand_name}_${solvent_name}.gro"
                echo "      Running: $solvate_cmd"
                solvate_out=$($solvate_cmd 2>&1 | tee solvate.log)
                # Parse the number of solvent molecules that were added
                solvent_mol=$(grep -oP 'Added\s+\K\d+(?=\s+molecules)' solvate.log)
                if [[ -z "$solvent_mol" ]]; then
                    echo "      WARNING: Could not parse solvent molecule count. Defaulting to 900."
                    solvent_mol=900
                fi
            fi
            echo "      Solvent molecules added: $solvent_mol"
            
            # Step 4: Create working topology files (both standard and 0q versions) from the original backup
            cp "../../${ligand_name}_GMX_original.top" "${ligand_name}_GMX.top"
            cp "${ligand_name}_GMX.top" "${ligand_name}_GMX_0q.top"
            echo "      Created working topology files from original backup."

            # For non-water solvents, we insert our own forcefield.itp lines and remove default headers
            # This is necessary to handle partial charges vs. zero-charges in separate .itp files
            if [[ "$solvent_name" != "H2O" && "$solvent_name" != "D2O" ]]; then
                rm -f ligand_atomtypes.tmp
                for top_file in ${ligand_name}_GMX.top ${ligand_name}_GMX_0q.top; do
                    atom_line=$(grep -n "^\[ *atomtypes *\]" "$top_file" | cut -d: -f1 | head -n1)
                    next_section=$(grep -n "^\[" "$top_file" | awk -F: -v start="$atom_line" '$1 > start {print $1; exit}')
                    if [[ -z "$next_section" ]]; then
                        next_section=$(wc -l < "$top_file")
                    fi
                    # Extract lines from [ atomtypes ] to the next section
                    sed -n "$((atom_line+1)),$((next_section - 1))p" "$top_file" >> ligand_atomtypes.tmp
                    
                    # Remove everything up to next_section - 1
                    sed -i "1,$((next_section - 1))d" "$top_file"
                    
                    # Now prepend new header block with includes
                    {
                      echo '#include "amber99sb.ff/forcefield.itp"'
                      echo '#include "./forcefield.itp"'
                      echo ""
                      cat "$top_file"
                    } > tmp_top && mv tmp_top "$top_file"
                done
                
                # Combine the ligand's atomtypes with any from topol.top to form forcefield.itp
                {
                  echo "[ atomtypes ]"
                  # Remove duplicates via awk + a small hash-based filter
                  awk '!/^ *$/ { if (!seen[$0]++) print }' ligand_atomtypes.tmp
                  if [[ -f topol.top ]]; then
                      solvent_atom_line=$(grep -n "^\[ *atomtypes *\]" topol.top | cut -d: -f1 | head -n1)
                      solvent_next_section=$(grep -n "^\[" topol.top | awk -F: -v start="$solvent_atom_line" '$1 > start {print $1; exit}')
                      if [[ -z "$solvent_next_section" ]]; then
                          solvent_next_section=$(wc -l < topol.top)
                      fi
                      sed -n "$((solvent_atom_line+1)),$((solvent_next_section - 1))p" topol.top | awk '!/^ *$/ { if (!seen[$0]++) print }'
                  fi
                } | awk '!/^ *$/ { if (!seen[$0]++) print }' > forcefield.itp
                rm -f ligand_atomtypes.tmp
            else
                # For H2O/D2O, we remove any header lines above [ atomtypes ] and prepend a simple header
                for top_file in ${ligand_name}_GMX.top ${ligand_name}_GMX_0q.top; do
                    atom_line=$(grep -n "^\[ *atomtypes *\]" "$top_file" | cut -d: -f1 | head -n1)
                    if [[ -n "$atom_line" && "$atom_line" -gt 1 ]]; then
                        sed -i "1,$((atom_line - 1))d" "$top_file"
                    fi
                    {
                      echo '#include "amber99sb.ff/forcefield.itp"'
                      echo ""
                      cat "$top_file"
                    } > tmp_top && mv tmp_top "$top_file"
                done
            fi
            
            # Step 5: Insert the solvent include directives above the [ system ] section
            # For water, we reference tip3p.itp or tip3p_0q.itp; for others we reference <solvent>.itp or <solvent>_0q.itp
            if [[ "$solvent_name" != "H2O" && "$solvent_name" != "D2O" ]]; then
                sed -i '/^\[ *system *\]/i\
\
#include ".\/'"$solvent_name"'.itp"\
\
' ${ligand_name}_GMX.top
                
                sed -i '/^\[ *system *\]/i\
\
#include ".\/'"$solvent_name"'_0q.itp"\
\
' ${ligand_name}_GMX_0q.top
            else
                # For H2O or D2O, we use the modified tip3p includes
                sed -i '/^\[ *system *\]/i\
\
#include ".\/tip3p.itp"\
\
' ${ligand_name}_GMX.top

                sed -i '/^\[ *system *\]/i\
\
#include ".\/tip3p_0q.itp"\
\
' ${ligand_name}_GMX_0q.top
            fi
            echo "      Modified topology headers and inserted include directives."

            # Step 6: Create solvent .itp files for non-water solvents
            # This means isolating the relevant lines from topol.top
            if [[ "$solvent_name" != "H2O" && "$solvent_name" != "D2O" ]]; then
                cp topol.top "${solvent_name}.itp"
                cp topol.top "${solvent_name}_0q.itp"

                # Retain only what's below the "[ moleculetype ]" section
                sed -i '1,/^\[ *moleculetype *\]/{/^\[ *moleculetype *\]/!d}' "${solvent_name}.itp"
                sed -i '1,/^\[ *moleculetype *\]/{/^\[ *moleculetype *\]/!d}' "${solvent_name}_0q.itp"

                # Remove [ system ] and [ molecules ] sections
                sed -i '/^\[ *system *\]/,/^$/d' "${solvent_name}.itp"
                sed -i '/^\[ *molecules *\]/,/^$/d' "${solvent_name}.itp"
                sed -i '/^\[ *system *\]/,/^$/d' "${solvent_name}_0q.itp"
                sed -i '/^\[ *molecules *\]/,/^$/d' "${solvent_name}_0q.itp"
                echo "      Created and cleaned solvent itp files."

                # Set charges to 0.0 in the [ atoms ] section for the _0q file
                awk 'BEGIN { in_atoms=0 }
                     /^\[ *atoms *\]/ { in_atoms=1; print; next }
                     in_atoms && /^\[/ { in_atoms=0 }
                     in_atoms && NF > 0 {
                         $7="0.0";
                         print; 
                         next
                     }
                     { print }' "${solvent_name}_0q.itp" > tmp_0q.itp && mv tmp_0q.itp "${solvent_name}_0q.itp"
                echo "      Set charges to 0.0 in ${solvent_name}_0q.itp"
            fi
            
            # Step 7: Rebuild the [ molecules ] section in both topology files
            # Remove any existing [ molecules ] block and add back with correct mol counts
            for top_file in ${ligand_name}_GMX.top ${ligand_name}_GMX_0q.top; do
                sed -i '/^\[ *molecules *\]/,$d' "$top_file"
                echo "" >> "$top_file"
                echo "[ molecules ]" >> "$top_file"
                echo "; Compound        nmols" >> "$top_file"
                echo " $ligand_name              1     " >> "$top_file"
                if [[ "$solvent_name" == "H2O" || "$solvent_name" == "D2O" ]]; then
                    echo "SOL              $solvent_mol" >> "$top_file"
                else
                    echo "MOL              $solvent_mol" >> "$top_file"
                fi
            done
            echo "      Created clean [ molecules ] section with correct solvent identifier."

            # Remove the temporary topology file now that the system is fully updated
            rm -f temp_topol.top

            # Step 8: Copy MD parameter files (min, nvt, npt, md) from MD_PARAMS
            cp "$MD_PARAMS_DIR/min.mdp" .
            cp "$MD_PARAMS_DIR/md.mdp" .
            cp "$MD_PARAMS_DIR/npt.mdp" .
            cp "$MD_PARAMS_DIR/nvt.mdp" .
            echo "      Copied mdp files from MD_PARAMS."

            # Step 9: Modify the mdp files for correct temperature coupling groups
            # We use the first 4 characters of the ligand name for the first T-coupling group
            ligand_tc_name="${ligand_name:0:4}"
            if [[ "$solvent_name" != "H2O" && "$solvent_name" != "D2O" ]]; then
                # For non-water, T-coupling groups are: <ligand_tc_name> MOL
                sed -i 's/tc-grps\s*=\s*#\s*#/tc-grps = '"${ligand_tc_name}"' MOL/' min.mdp md.mdp npt.mdp nvt.mdp
            else
                # For water or D2O, T-coupling groups are: <ligand_tc_name> SOL
                sed -i 's/tc-grps\s*=\s*#\s*#/tc-grps = '"${ligand_tc_name}"' SOL/' min.mdp md.mdp npt.mdp nvt.mdp
            fi
            echo "      Modified mdp files with correct tc-grps settings (using 4-char limit for ligand name)."

            # Copy all analysis scripts (if any) into the run folder
            if [[ -d "$ANALYSIS_SCRIPTS_DIR" ]]; then
                file_count=$(find "$ANALYSIS_SCRIPTS_DIR" -type f | wc -l)
                if [[ $file_count -gt 0 ]]; then
                    cp "$ANALYSIS_SCRIPTS_DIR"/* .
                    echo "      Copied $file_count analysis scripts from ANALYSIS_SCRIPTS directory"
                else
                    echo "      WARNING: ANALYSIS_SCRIPTS directory exists but contains no files"
                fi
            fi

            # Step 10: Create a submission script (run_md.sh) for MD
            # This script includes minimization, NVT, NPT, production, and zero-charge (0q) rerun
            cat > run_md.sh <<EOF
#!/bin/bash
## SLURM submission headers
#SBATCH --job-name=${solvent_name}_run${run_num}
#SBATCH --ntasks-per-node=4             # request 24 MPI tasks per node
#SBATCH --cpus-per-task=2                # 2 OpenMP threads per MPI task => total: 24 x 2 = 48 CPUs/node
#SBATCH --constraint="[skylake|cascade]" # restrict to AVX512 capable nodes.
#SBATCH --mem-per-cpu=2000M              # memory per CPU (in MB)
#SBATCH --time=0-01:00                   # time limit (D-HH:MM)
module purge
module load arch/avx512   # switch architecture for up to 30% speedup
module load  StdEnv/2023  gcc/12.3  openmpi/4.1.5  gromacs/2024.4
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Energy minimization using min.mdp
gmx grompp -f min.mdp -c ${ligand_name}_${solvent_name}.gro -p ${ligand_name}_GMX.top -o ${ligand_name}_${solvent_name}_min.tpr -maxwarn 2
gmx mdrun -v -deffnm ${ligand_name}_${solvent_name}_min

# NVT equilibration
gmx grompp -f nvt.mdp -c ${ligand_name}_${solvent_name}_min.gro -p ${ligand_name}_GMX.top -o ${ligand_name}_${solvent_name}_nvt.tpr -maxwarn 2
gmx mdrun -v -deffnm ${ligand_name}_${solvent_name}_nvt

# NPT equilibration
gmx grompp -f npt.mdp -c ${ligand_name}_${solvent_name}_nvt.gro -p ${ligand_name}_GMX.top -o ${ligand_name}_${solvent_name}_npt.tpr -maxwarn 2
gmx mdrun -v -deffnm ${ligand_name}_${solvent_name}_npt

# Production MD run
gmx grompp -f md.mdp -c ${ligand_name}_${solvent_name}_npt.gro -p ${ligand_name}_GMX.top -o ${ligand_name}_${solvent_name}_md.tpr -maxwarn 2
gmx mdrun -v -deffnm ${ligand_name}_${solvent_name}_md

# Trajectory analysis
gmx traj -s ${ligand_name}_${solvent_name}_md.tpr -f ${ligand_name}_${solvent_name}_md.trr -n ../../probe.ndx -ox co_coords.xvg
gmx traj -s ${ligand_name}_${solvent_name}_md.tpr -f ${ligand_name}_${solvent_name}_md.trr -n ../../probe.ndx -of co_forces.xvg

# Zero charge (0q) topology run
gmx grompp -f md.mdp -c ${ligand_name}_${solvent_name}_npt.gro -p ${ligand_name}_GMX_0q.top -o ${ligand_name}_${solvent_name}_md_0q.tpr
gmx mdrun -rerun ${ligand_name}_${solvent_name}_md.trr -v -deffnm ${ligand_name}_${solvent_name}_md_0q

gmx traj -s ${ligand_name}_${solvent_name}_md_0q.tpr -f ${ligand_name}_${solvent_name}_md_0q.trr -n ../../probe.ndx -ox co_coords_0q.xvg
gmx traj -s ${ligand_name}_${solvent_name}_md_0q.tpr -f ${ligand_name}_${solvent_name}_md_0q.trr -n ../../probe.ndx -of co_forces_0q.xvg
EOF

            chmod +x run_md.sh
            echo "      Created MD submission script: run_md.sh with job name ${solvent_name}_run${run_num}"
            
            # Return to the ligand directory for the next run
            cd "$ligand_abs" || continue
        done
        
        echo "    Completed triplicate setup for solvent: $solvent_name in ligand: $ligand_name"
    done
    
    # After all solvents are processed, restore the original topology
    if [[ -f "${ligand_name}_GMX_original.top" ]]; then
        mv "${ligand_name}_GMX_original.top" "${ligand_name}_GMX.top"
        echo "Restored original topology file for ligand: $ligand_name"
    fi
    
    echo "Completed MD setup for ligand: $ligand_name"
    cd "$script_dir" || exit
done

echo "=============================================="
echo "Triplicate MD setup for all ligands and solvents has been completed."
echo "Analysis scripts have been copied to each run directory."
