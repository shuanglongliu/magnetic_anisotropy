import os

root_dir = os.path.dirname(os.path.abspath(__file__))  + "/"

incar = """
SYSTEM = vasp

#### symmetry ####
#ISYM = 0

#### system size ####
#NBANDS = 448
#NELECT = 664.0

#### accuracy ####
PREC = Accurate
ENCUT = 600
LREAL = F

#### parallelization ####
#KPAR = 4
NCORE = 16

#### density functional ####

#=== LDA/GGA ===#
#GGA = PE

#=== METAGGA ===#
#METAGGA = R2SCAN

#=== Hybrid: HSE06 ===#
#LHFCALC = T 
#GGA = PE
#HFSCREEN = 0.2 
#PRECFOCK = Accurate

#=== vdW-DF: optB86b ===#
#GGA      = MK 
#PARAM1   = 0.1234 
#PARAM2   = 1.0
#AGGAC    = 0.0
#LUSE_VDW = .TRUE.
#LASPH    = .TRUE.

#=== LIBXC ===#
#GGA = LIBXC
#LIBXC1 = GGA_X_PBE
#LIBXC2 = GGA_C_PBE

#### empirical vdW ####
#IVDW = 11

#### LDA+U ####
#LDAU = T
#LDAUTYPE = 2
#LDAUPRINT = 1
#LDAUL =   2   -1    -1   -1   -1   -1   -1  
#LDAUU =   2.5  0.0   0.0  0.0  0.0  0.0  0.0
#LDAUJ =   0.0  0.0   0.0  0.0  0.0  0.0  0.0

#### electronic optimization ####
EDIFF = 1E-8
NELM = 120
#NELMIN = 10
#NELMDL = -12
AMIX = 0.2
AMIX_MAG = 0.4
AMIX_MIN = 0.05
BMIX = 0.0001
BMIX_MAG = 0.0001
ALGO = All
TIME = 0.35

#### structural relaxation ####
#NSW = 500
#IBRION = 2
#ISIF = 2
#EDIFFG = -0.02
#POTIM = 0.3

#### magnetism: accuracy ####
LASPH = T
GGA_COMPAT = F

#### magnetism: collinear spin ####
#ISPIN = 2
#MAGMOM = 3.0 3.0 3.0 1000*0.0
#NUPDOWN = 9

#### magnetism: noncollinear spin, SOC ####
LSORBIT = T
SAXIS = 0 0 1
MAGMOM = {mstring:s} 3000*0.0

#### magnetism: constraint ####
I_CONSTRAINED_M = 1
M_CONSTR = {mstring:s} 3000*0.0

LAMBDA =  15.0
RWIGS = {rwigs:s}

#### magnetism: orbital moment ####
LORBMOM = T

#### charge, wavefunction ####
ISTART = 1
ICHARG = 0
LWAVE = T
LCHARG = F
LAECHG = F
LMAXMIX = 4

#### dos ####
ISMEAR = 0
SIGMA = 0.01
#NEDOS = 2501
#EMIN = -15
#EMAX = 10
LORBIT = 10

#### Partial Charge ####
#LPARD = T
#EINT = -5.0000  -4.7355
#EINT = -4.7355 -4.0500

#### wann ####
#LWANNIER90 = .T.
#LWRITE_UNK = .TRUE.

#### polarization ####
#IDIPOL = 4
#LDIPOL = T
#DIPOL = 0.5 0.5 0.5
#LCALCPOL = T

#### occupation matrix control ####
#OCCEXT = 1
"""

kpoints = """kpoints
0
gamma
 1   1   1
 0.0 0.0 0.0
"""

poscar = """
"""

job_script_knl_mpi = """#!/bin/bash

#SBATCH -A m3346
#SBATCH -J Mn4Na
#SBATCH -q regular
#SBATCH -C knl 
#SBATCH -N 5
#SBATCH --ntasks-per-node=64
#SBATCH --time=0-6:00:00 
#SBATCH --error=error
#SBATCH --output=output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shlufl@ufl.edu
##SBATCH --dependency=afterok:45514509
 
module load vasp/6.3.2-knl

export OMP_NUM_THREADS=1

keep_log

srun -n $SLURM_NTASKS -c 4 --cpu_bind=cores vasp_ncl 

rm CHG CHGCAR PROCAR DOSCAR vasprun.xml vaspout.h5
sed -i "/PROFILE/d" output
"""

job_script_knl_hybrid = """#!/bin/bash

#SBATCH -A m3346
#SBATCH -J Co3
#SBATCH -q regular
#SBATCH -C knl
#SBATCH -N 4
#SBATCH --ntasks-per-node=8
#SBATCH --time=1-00:00:00
#SBATCH --error=error
#SBATCH --output=output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shlufl@ufl.edu
##SBATCH --dependency=afterok:45637525

module load vasp/6.1.2-knl

export OMP_NUM_THREADS=8

keep_log

srun -n 32 -c 32 --cpu_bind=cores vasp_ncl
"""

job_script_hipergator_rome = """#!/bin/bash -l

#SBATCH --account=m2qm-efrc
#SBATCH --qos=m2qm-efrc-b
#SBATCH --job-name=Co3
#SBATCH --mail-type=All
#SBATCH --mail-user=shlufl@ufl.edu
#SBATCH --partition=hpg-default
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-socket=16
#SBATCH --mem=200gb
#SBATCH --distribution=cyclic:cyclic
#SBATCH -t 96:00:00
#SBATCH --err=error
#SBATCH --output=output
##SBATCH --dependency=afterok:8190473

module purge
module load intel/2020.0.166 openmpi/4.1.1

source ~/.bash_aliases; keep_log

srun --mpi=pmix_v2  $HOME/bin/vasp_ncl
"""

job_script_hipergator_a100 = """#!/bin/bash

#SBATCH --job-name=Cr3
#SBATCH --mem-per-cpu=32gb
#SBATCH -t 4-00:00:00
#SBATCH -p gpu --gpus=a100:2
#SBATCH --account=m2qm-efrc
#SBATCH --qos=m2qm-efrc
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shlufl@ufl.edu
#SBATCH --error=error
#SBATCH --output=output
##SBATCH --dependency=afterok:9106099

module purge; module load cuda/11.1.0 nvhpc/20.11 openmpi/4.0.5 qd/2.3.22 fftw/3.3.8 vasp/6.2.0

source ~/.bash_aliases; keep_log

srun --mpi=pmix vasp_ncl
"""

job_script_perlmutter_cpu = """#!/bin/bash

#SBATCH -A m3346
#SBATCH -J Mn4
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -n 128
#SBATCH --error=error
#SBATCH --output=output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shlufl@ufl.edu
##SBATCH --dependency=afterok:3760298

module load PrgEnv-nvidia/8.3.3 vasp/6.3.2-cpu

export OMP_NUM_THREADS=1

source ~/.bash_aliases; keep_log

srun -n 128 -c 2 --cpu-bind=cores vasp_ncl

rm -rf CHG CHGCAR DOSCAR PROCAR vasprun.xml vaspout.h5

sed -i "/PROFILE/d" output
"""

job_script_perlmutter_gpu_1_node = """#!/bin/bash

#SBATCH -A m3346_g
#SBATCH -J Mn4Na
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 32
#SBATCH -G 4
#SBATCH --exclusive
#SBATCH --error=error
#SBATCH --output=output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shlufl@ufl.edu
##SBATCH --dependency=afterok:45514509

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

module load vasp/6.3.2-gpu

source ~/.bash_aliases; keep_log

srun -n 4 -c 32 --cpu-bind=cores --gpu-bind=none vasp_ncl

rm -rf CHG CHGCAR DOSCAR PROCAR vasprun.xml vaspout.h5

sed -i "/PROFILE/d" output
"""

job_script_perlmutter_gpu_4_nodes = """#!/bin/bash

#SBATCH -A m3346_g
#SBATCH -J FePc
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 02:00:00
#SBATCH -N 4
#SBATCH -n 16
#SBATCH -c 32
#SBATCH -G 16
#SBATCH --exclusive
#SBATCH --error=error
#SBATCH --output=output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shlufl@ufl.edu
##SBATCH --dependency=afterok:45514509

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

module load vasp/6.3.2-gpu

source ~/.bash_aliases; keep_log

srun -n 16 -c 32 --cpu-bind=cores --gpu-bind=none vasp_ncl

rm -rf CHG vasprun.xml vaspout.h5

sed -i "/PROFILE/d" output
"""

job_script = job_script_perlmutter_gpu_1_node

emats_file = []

for i in range(100):
    emats_file.append( [[1,0,0],[0,1,0],[0,0,1]] )

