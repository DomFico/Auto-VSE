#!/bin/bash
## SLURM submission headers
#SBATCH --job-name=trichloromethane_run1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=0:30:00
module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 gromacs/2024.4
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Energy minimization using min.mdp
gmx grompp -f min.mdp -c METB_trichloromethane.gro -p METB_GMX.top -o METB_trichloromethane_min.tpr -maxwarn 2
gmx mdrun -v -deffnm METB_trichloromethane_min

# NVT equilibration
gmx grompp -f nvt.mdp -c METB_trichloromethane_min.gro -p METB_GMX.top -o METB_trichloromethane_nvt.tpr -maxwarn 2
gmx mdrun -v -deffnm METB_trichloromethane_nvt -nt 24

# NPT equilibration
gmx grompp -f npt.mdp -c METB_trichloromethane_nvt.gro -p METB_GMX.top -o METB_trichloromethane_npt.tpr -maxwarn 2
gmx mdrun -v -deffnm METB_trichloromethane_npt -nt 24

# Production MD run
gmx grompp -f md.mdp -c METB_trichloromethane_npt.gro -p METB_GMX.top -o METB_trichloromethane_md.tpr -maxwarn 2
gmx mdrun -v -deffnm METB_trichloromethane_md -nt 24

# Trajectory analysis
gmx traj -s METB_trichloromethane_md.tpr -f METB_trichloromethane_md.trr -n ../../probe.ndx -ox co_coords.xvg
gmx traj -s METB_trichloromethane_md.tpr -f METB_trichloromethane_md.trr -n ../../probe.ndx -of co_forces.xvg

# Zero charge (0q) topology run
gmx grompp -f md.mdp -c METB_trichloromethane_npt.gro -p METB_GMX_0q.top -o METB_trichloromethane_md_0q.tpr
gmx mdrun -rerun METB_trichloromethane_md.trr -v -deffnm METB_trichloromethane_md_0q

gmx traj -s METB_trichloromethane_md_0q.tpr -f METB_trichloromethane_md_0q.trr -n ../../probe.ndx -ox co_coords_0q.xvg
gmx traj -s METB_trichloromethane_md_0q.tpr -f METB_trichloromethane_md_0q.trr -n ../../probe.ndx -of co_forces_0q.xvg
