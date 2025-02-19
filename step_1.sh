#!/bin/bash
# -----------------------------------------------------------------------------
# step_1.sh
#
# This script processes each ligand subdirectory in the LIGANDS folder.
# It reads the SMILES string from a file named smiles.txt, converts that
# to a mol2 file, assigns GAFF atom types and force field parameters,
# runs tleap to generate AMBER topology files, converts them to GROMACS
# format with acepype, and organizes all outputs.
#
# IMPORTANT NOTE: This protocol only deals with neutral compounds.
#                 If a "step_1_outputs" folder exists, the user can
#                 decide to re-run by cleaning old outputs.
#
# Steps (for each ligand folder):
#   1) Check for and remove any unwanted directories (*_md or "#*").
#   2) Confirm if step_1_outputs already exists; if so, prompt user to re-run.
#   3) Parse the SMILES from smiles.txt.
#   4) Convert SMILES to mol2 (using obabel).
#   5) Assign GAFF atom types (antechamber).
#   6) Generate force field parameters (parmchk2).
#   7) Produce AMBER topologies (tleap).
#   8) Convert AMBER files to GROMACS format (acepype).
#   9) Move all intermediates into step_1_outputs folder and keep final .gro/.top.
# -----------------------------------------------------------------------------

# Iterate through each subdirectory under LIGANDS/
for ligand_dir in LIGANDS/*/; do
    
    # Extract directory name without trailing slash
    ligand_name=$(basename "$ligand_dir")
    
    echo "=============================================="
    echo "Processing ligand: $ligand_name"

    # Enter the ligand directory or skip if unable
    cd "$ligand_dir" || { echo "ERROR: Cannot cd into $ligand_dir"; continue; }
    
    # -------------------------------------------------------------------------
    # Remove any unwanted directories that match *_md or "#*"
    # This helps keep the folder clean if leftover MD runs exist.
    # -------------------------------------------------------------------------
    for bad_dir in *_md "#"*; do
        if [[ -d "$bad_dir" ]]; then
            echo "Removing unwanted directory: $bad_dir"
            rm -rf "$bad_dir"
        fi
    done

    # -------------------------------------------------------------------------
    # If 'step_1_outputs' folder already exists, ask the user if they want to 
    # re-run. If yes, we remove prior outputs; if no, skip this ligand.
    # -------------------------------------------------------------------------
    if [[ -d step_1_outputs ]]; then
        read -p "Ligand '$ligand_name': A step_1_outputs folder exists. Re-run the process? [y/n]: " run_again
        if [[ "$run_again" =~ ^[Yy]$ ]]; then
            echo "Cleaning existing step_1_outputs and intermediate files for $ligand_name..."
            rm -rf step_1_outputs \
                   "${ligand_name}.mol2" \
                   "${ligand_name}_gaff.mol2" \
                   "${ligand_name}.frcmod" \
                   "${ligand_name}.prmtop" \
                   "${ligand_name}.inpcrd" \
                   tleap.in \
                   "${ligand_name}.lib" \
                   "${ligand_name}.amb2gmx"
        else
            echo "Skipping ligand: $ligand_name"
            cd - > /dev/null
            continue
        fi
    fi

    # -------------------------------------------------------------------------
    # Check if smiles.txt exists. If not present, we cannot proceed.
    # -------------------------------------------------------------------------
    if [[ ! -f "smiles.txt" ]]; then
        echo "ERROR: smiles.txt not found in $ligand_name. Skipping this ligand."
        cd - > /dev/null
        continue
    fi

    # -------------------------------------------------------------------------
    # Create the step_1_outputs folder. This will hold logs and intermediate files.
    # -------------------------------------------------------------------------
    mkdir -p step_1_outputs

    # Initialize the command log file
    CMD_LOG="step_1_outputs/commands.txt"
    echo "Commands used for ligand $ligand_name:" > "$CMD_LOG"

    # =========================================================================
    # Step 1: Parse SMILES from the first line of smiles.txt
    # =========================================================================
    SMILES=$(head -n 1 smiles.txt)
    echo "Parsed SMILES: $SMILES"
    echo "SMILES: $SMILES" >> "$CMD_LOG"

    # =========================================================================
    # Step 2: Convert the SMILES to a mol2 file using obabel with 3D geometry
    # =========================================================================
    OBA_CMD="obabel -:\"$SMILES\" -o mol2 -O ${ligand_name}.mol2 --gen3d --best"
    echo "Running: $OBA_CMD"
    echo "$OBA_CMD" >> "$CMD_LOG"
    eval $OBA_CMD

    # =========================================================================
    # Step 3: Use antechamber to assign GAFF atom types
    # Note: This step assumes neutral molecules only
    # =========================================================================
    AC_CMD="antechamber -i ${ligand_name}.mol2 -fi mol2 -o ${ligand_name}_gaff.mol2 -fo mol2 -c bcc -s 2"
    echo "Running: $AC_CMD"
    echo "$AC_CMD" >> "$CMD_LOG"
    eval $AC_CMD

    # =========================================================================
    # Step 4: Use parmchk2 to produce the GAFF force field parameters (.frcmod)
    # =========================================================================
    PCMD="parmchk2 -i ${ligand_name}_gaff.mol2 -f mol2 -o ${ligand_name}.frcmod"
    echo "Running: $PCMD"
    echo "$PCMD" >> "$CMD_LOG"
    eval $PCMD

    # =========================================================================
    # Step 5: Generate AMBER topology files with tleap
    #   - We create a custom tleap.in file, then run it
    # =========================================================================
    cat > tleap.in <<EOF
source leaprc.gaff
${ligand_name} = loadmol2 ${ligand_name}_gaff.mol2
loadamberparams ${ligand_name}.frcmod
check ${ligand_name}
saveoff ${ligand_name} ${ligand_name}.lib
saveamberparm ${ligand_name} ${ligand_name}.prmtop ${ligand_name}.inpcrd
quit
EOF
    echo "Created tleap.in file." >> "$CMD_LOG"
    
    TLEAP_CMD="tleap -f tleap.in"
    echo "Running: $TLEAP_CMD"
    echo "$TLEAP_CMD" >> "$CMD_LOG"
    eval $TLEAP_CMD

    # =========================================================================
    # Step 6: Modify the .prmtop file to replace 'UNL' with the ligand name
    # =========================================================================
    SED_CMD="sed -i 's/UNL/${ligand_name}/g' ${ligand_name}.prmtop"
    echo "Running: $SED_CMD"
    echo "$SED_CMD" >> "$CMD_LOG"
    eval $SED_CMD

    # =========================================================================
    # Step 6 (continued): Convert AMBER files to GROMACS format using acepype
    # =========================================================================
    ACPYPE_CMD="acpype -p ${ligand_name}.prmtop -x ${ligand_name}.inpcrd"
    echo "Running: $ACPYPE_CMD"
    echo "$ACPYPE_CMD" >> "$CMD_LOG"
    eval $ACPYPE_CMD

    # =========================================================================
    # Step 7: Retrieve final GROMACS files from the acepype output folder.
    #         The folder is named ${ligand_name}.amb2gmx by default.
    # =========================================================================
    AMB2GMX_DIR="${ligand_name}.amb2gmx"
    if [[ -d "$AMB2GMX_DIR" ]]; then
        MV_GRO="mv ${AMB2GMX_DIR}/${ligand_name}_GMX.gro ."
        MV_TOP="mv ${AMB2GMX_DIR}/${ligand_name}_GMX.top ."
        echo "Running: $MV_GRO"
        echo "$MV_GRO" >> "$CMD_LOG"
        eval $MV_GRO
        echo "Running: $MV_TOP"
        echo "$MV_TOP" >> "$CMD_LOG"
        eval $MV_TOP
    else
        echo "WARNING: Acepype output directory $AMB2GMX_DIR not found." | tee -a "$CMD_LOG"
    fi

    # =========================================================================
    # Step 7 (continued): Move all intermediate files to the step_1_outputs folder
    # leaving the final .gro and .top in the current directory.
    # Also keep smiles.txt and probe.ndx for future reprocessing.
    # =========================================================================
    echo "Organizing intermediate files into step_1_outputs folder."
    for f in *; do
        if [[ "$f" != *"_GMX.gro" && "$f" != *"_GMX.top" && "$f" != "step_1_outputs" && "$f" != "smiles.txt" && "$f" != "probe.ndx" ]]; then
            mv "$f" step_1_outputs/
        fi
    done
    # Ensure final .gro/.top remain at the top level if they got moved.
    if [[ -f "step_1_outputs/${ligand_name}_GMX.gro" ]]; then
        mv "step_1_outputs/${ligand_name}_GMX.gro" .
    fi
    if [[ -f "step_1_outputs/${ligand_name}_GMX.top" ]]; then
        mv "step_1_outputs/${ligand_name}_GMX.top" .
    fi

    # =========================================================================
    # Append a detailed, numbered list of commands with their functions to
    # the command log (commands.txt)
    # =========================================================================
    cat >> "step_1_outputs/commands.txt" <<EOF

Detailed Command List:
1. Parse SMILES:
   Command: head -n 1 smiles.txt
   Function: Extract the SMILES string from the smiles.txt file.

2. Convert SMILES to mol2:
   Command: obabel -:"$SMILES" -o mol2 -O ${ligand_name}.mol2 --gen3d --best
   Function: Generate a mol2 file with optimal geometry using a systematic rotamer search.

3. Assign GAFF Atom Types:
   Command: antechamber -i ${ligand_name}.mol2 -fi mol2 -o ${ligand_name}_gaff.mol2 -fo mol2 -c bcc -s 2
   Function: Assign GAFF atom types to the molecule (currently supports only neutral compounds).

4. Generate Force Field Parameters:
   Command: parmchk2 -i ${ligand_name}_gaff.mol2 -f mol2 -o ${ligand_name}.frcmod
   Function: Create a force field parameter file (frcmod) using GAFF parameters.

5. Generate Topology Files using tleap:
   Command: tleap -f tleap.in
   Function: Run tleap with a custom input file to generate AMBER topology files (.prmtop, .inpcrd, .lib).

6. Modify Residue Name in Topology:
   Command: sed -i 's/UNL/${ligand_name}/g' ${ligand_name}.prmtop
   Function: Replace the default residue name "UNL" with the ligand's directory name in the prmtop file.

7. Convert AMBER Files to GROMACS Format:
   Command: acpype -p ${ligand_name}.prmtop -x ${ligand_name}.inpcrd
   Function: Convert AMBER topology and coordinate files to GROMACS format (.gro, .top).

8. Retrieve Final Files:
   Command: mv ${ligand_name}.amb2gmx/${ligand_name}_GMX.gro .  and  mv ${ligand_name}.amb2gmx/${ligand_name}_GMX.top .
   Function: Extract the final GROMACS .gro and .top files from the acepype output directory.

9. Organize Intermediate Files:
   Command: Move all non-final files (except smiles.txt and probe.ndx) into the "step_1_outputs" folder.
   Function: Clean up the working directory by relocating intermediate files while preserving smiles.txt and probe.ndx.
EOF

    echo "Finished processing ligand: $ligand_name"

    # Return to the parent directory (the main project folder)
    cd - > /dev/null
done

echo "=============================================="
echo "All ligand directories have been processed."
