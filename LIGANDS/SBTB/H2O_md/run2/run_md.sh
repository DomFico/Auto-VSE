#!/bin/bash
## SLURM submission headers
#SBATCH --job-name=H2O_run2
#SBATCH --ntasks-per-node=4             # request 24 MPI tasks per node
#SBATCH --cpus-per-task=2                # 2 OpenMP threads per MPI task => total: 24 x 2 = 48 CPUs/node
#SBATCH --constraint="[skylake|cascade]" # restrict to AVX512 capable nodes.
#SBATCH --mem-per-cpu=2000M              # memory per CPU (in MB)
#SBATCH --time=0-01:00                   # time limit (D-HH:MM)
module purge
module load arch/avx512   # switch architecture for up to 30% speedup
module load  StdEnv/2023  gcc/12.3  openmpi/4.1.5  gromacs/2024.4
export OMP_NUM_THREADS="1"

# Energy minimization using min.mdp
gmx grompp -f min.mdp -c SBTB_H2O.gro -p SBTB_GMX.top -o SBTB_H2O_min.tpr -maxwarn 2
gmx mdrun -v -deffnm SBTB_H2O_min

# NVT equilibration
gmx grompp -f nvt.mdp -c SBTB_H2O_min.gro -p SBTB_GMX.top -o SBTB_H2O_nvt.tpr -maxwarn 2
gmx mdrun -v -deffnm SBTB_H2O_nvt

# NPT equilibration
gmx grompp -f npt.mdp -c SBTB_H2O_nvt.gro -p SBTB_GMX.top -o SBTB_H2O_npt.tpr -maxwarn 2
gmx mdrun -v -deffnm SBTB_H2O_npt

# Production MD run
gmx grompp -f md.mdp -c SBTB_H2O_npt.gro -p SBTB_GMX.top -o SBTB_H2O_md.tpr -maxwarn 2
gmx mdrun -v -deffnm SBTB_H2O_md

# Trajectory analysis
gmx traj -s SBTB_H2O_md.tpr -f SBTB_H2O_md.trr -n ../../probe.ndx -ox co_coords.xvg
gmx traj -s SBTB_H2O_md.tpr -f SBTB_H2O_md.trr -n ../../probe.ndx -of co_forces.xvg

# Zero charge (0q) topology run
gmx grompp -f md.mdp -c SBTB_H2O_npt.gro -p SBTB_GMX_0q.top -o SBTB_H2O_md_0q.tpr
gmx mdrun -rerun SBTB_H2O_md.trr -v -deffnm SBTB_H2O_md_0q

gmx traj -s SBTB_H2O_md_0q.tpr -f SBTB_H2O_md_0q.trr -n ../../probe.ndx -ox co_coords_0q.xvg
gmx traj -s SBTB_H2O_md_0q.tpr -f SBTB_H2O_md_0q.trr -n ../../probe.ndx -of co_forces_0q.xvg
