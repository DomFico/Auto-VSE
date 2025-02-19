#!/bin/bash
## SLURM submission headers
#SBATCH --job-name=D2O_run3
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=0:30:00
module purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 gromacs/2024.4
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Energy minimization using min.mdp
gmx grompp -f min.mdp -c E4TPB_D2O.gro -p E4TPB_GMX.top -o E4TPB_D2O_min.tpr -maxwarn 2
gmx mdrun -v -deffnm E4TPB_D2O_min

# NVT equilibration
gmx grompp -f nvt.mdp -c E4TPB_D2O_min.gro -p E4TPB_GMX.top -o E4TPB_D2O_nvt.tpr -maxwarn 2
gmx mdrun -v -deffnm E4TPB_D2O_nvt -nt 24

# NPT equilibration
gmx grompp -f npt.mdp -c E4TPB_D2O_nvt.gro -p E4TPB_GMX.top -o E4TPB_D2O_npt.tpr -maxwarn 2
gmx mdrun -v -deffnm E4TPB_D2O_npt -nt 24

# Production MD run
gmx grompp -f md.mdp -c E4TPB_D2O_npt.gro -p E4TPB_GMX.top -o E4TPB_D2O_md.tpr -maxwarn 2
gmx mdrun -v -deffnm E4TPB_D2O_md -nt 24

# Trajectory analysis
gmx traj -s E4TPB_D2O_md.tpr -f E4TPB_D2O_md.trr -n ../../probe.ndx -ox co_coords.xvg
gmx traj -s E4TPB_D2O_md.tpr -f E4TPB_D2O_md.trr -n ../../probe.ndx -of co_forces.xvg

# Zero charge (0q) topology run
gmx grompp -f md.mdp -c E4TPB_D2O_npt.gro -p E4TPB_GMX_0q.top -o E4TPB_D2O_md_0q.tpr
gmx mdrun -rerun E4TPB_D2O_md.trr -v -deffnm E4TPB_D2O_md_0q

gmx traj -s E4TPB_D2O_md_0q.tpr -f E4TPB_D2O_md_0q.trr -n ../../probe.ndx -ox co_coords_0q.xvg
gmx traj -s E4TPB_D2O_md_0q.tpr -f E4TPB_D2O_md_0q.trr -n ../../probe.ndx -of co_forces_0q.xvg
