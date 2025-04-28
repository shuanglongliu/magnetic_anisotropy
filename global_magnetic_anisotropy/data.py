root_dir = "/global/cfs/cdirs/m3346/shlufl/IETS/nco_basecell_beta_strain/mae/"
wave_dir = "6476401"

incar = """
SYSTEM = vasp

#### sym ####
#ISYM = 0

#### system size ####
#NBANDS = 448
#NELECT = 383.0

#### accuracy ####
PREC = Accurate
ENCUT = 450
LREAL = F

#### parallelization ####
KPAR = 4
NCORE = 8

#### electronic optimization ####
EDIFF = 1E-8
NELM = 90
#NELMIN = 10
#NELMDL = -8
AMIX = 0.2
AMIX_MAG = 0.4
AMIX_MIN = 0.05
BMIX = 0.0001
BMIX_MAG = 0.0001
ALGO = All

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
SAXIS = {sstring:s}
MAGMOM = 0.0 0.0 -1.0 \\
         0.0 0.0 -1.0 \\
         0.0 0.0 -1.0 \\
         0.0 0.0 -1.0 \\
         0.0 0.0 3.0 \\
         0.0 0.0 3.0 \\
         0.0 0.0 3.0 \\
         0.0 0.0 3.0 \\
         3300*0.0

#### magnetism: constraint ####
I_CONSTRAINED_M = 1
M_CONSTR =  0.0 0.0 -1.0 \\
         0.0 0.0 -1.0 \\
         0.0 0.0 -1.0 \\
         0.0 0.0 -1.0 \\
         0.0 0.0 3.0 \\
         0.0 0.0 3.0 \\
         0.0 0.0 3.0 \\
         0.0 0.0 3.0 \\
         3300*0.0

LAMBDA =  10.0
RWIGS = 1.286 1.302 0.820

#### magnetism: orbital moment ####
LORBMOM = T

#### charge, wavefunction ####
ISTART = 1
ICHARG = 0
LWAVE = F
LCHARG = F
LAECHG = F
LMAXMIX = 4

#### dos ####
ISMEAR = 0
SIGMA = 0.02
NEDOS = 501
EMIN = -15
EMAX = 10
LORBIT = 11

#### Partial Charge ####
#LPARD = T
#EINT = -5.0000  -4.7355
#EINT = -4.7355 -4.0500

#### vdW ####
IVDW = 11

#### LDA+U ####
LDAU = T
LDAUTYPE = 2
LDAUPRINT = 1
LDAUL =   2    2   -1  
LDAUU =   5.5  3.0  0.0
LDAUJ =   0.0  0.0  0.0

#### HSE ####
#LHFCALC = T 
#HFSCREEN = 0.2 
#PRECFOCK = Accurate
#ALGO = All 
#TIME = 0.35

#### wann ####
#LWANNIER90 = .T.
#LWRITE_UNK = .TRUE.

### polarization ###
#EFIELD_PEAD = 0.000 0.000 0.000
#IDIPOL = 3
#LMONO = T
#LDIPOL = T
#LCALCPOL = T
#DIPOL = 0.5 0.5 0.5 

### occupation matrix control ###
#OCCEXT = 1
"""

kpoints = """kpoints
0
gamma
 10  10  8
 0.0 0.0 0.0
"""

poscar = """
"""

job_script_mpi = """#!/bin/bash

#SBATCH -A m3346
#SBATCH -J Co3
#SBATCH -q regular
#SBATCH -C knl 
#SBATCH -N 2
#SBATCH --ntasks-per-node=64
#SBATCH --time=0-10:00:00 
#SBATCH --error=error
#SBATCH --output=output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shlufl@ufl.edu
##SBATCH --dependency=afterok:45514509
 
module load vasp/6.2.1-knl

export OMP_NUM_THREADS=1

keep_log

srun -n 128 -c 4 --cpu_bind=cores vasp_ncl 
"""

job_script_hybrid = """#!/bin/bash

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

job_script_hipergator = """#!/bin/bash -l

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

job_script_hpg_gpu = """#!/bin/bash

#SBATCH --job-name=Cr3
#SBATCH --mem-per-cpu=20gb
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

module purge; module load cuda/11.1.0 nvhpc/20.11 qd/2.3.22 openmpi/4.0.5 fftw/3.3.8

export LD_LIBRARY_PATH=/apps/nvidia/nvhpc/Linux_x86_64/20.11/compilers/extras/qd/lib:$LD_LIBRARY_PATH

source ~/.bash_aliases; keep_log

srun --mpi=pmix /home/shlufl/apps/vasp.6.2.1/bin_gpu/vasp_ncl
"""

job_script_perlmutter_1_node = """#!/bin/bash
#SBATCH -A m3346_g
#SBATCH -J Cr3
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -G 4
#SBATCH -c 32
#SBATCH --gpus-per-task=1
#SBATCH --ntasks-per-node=4
#SBATCH --error=error
#SBATCH --output=output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shlufl@ufl.edu
##SBATCH --dependency=afterok:45514509

module load vasp/6.2.1-gpu

source ~/.bash_aliases; keep_log

#if [ -f CONTCAR ]; then nl=$(wc -l CONTCAR | awk '{printf "%d", $1}'); if [ $nl -gt 8 ]; then cp CONTCAR POSCAR; else echo "Incomplete CONTCAR"; exit; fi; fi
#if [ -f STOPCAR ]; then rm STOPCAR; fi

srun -n 4 -c 32 --cpu-bind=cores -G 4 --gpu-bind=single:1 vasp_ncl
"""

job_script_perlmutter_4_nodes = """#!/bin/bash

#SBATCH -A m3346_g
#SBATCH -J NCO
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 02:00:00
#SBATCH -N 4
#SBATCH -n 16
#SBATCH -G 16
#SBATCH -c 32
#SBATCH --gpus-per-task=1
#SBATCH --error=error
#SBATCH --output=output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shlufl@ufl.edu
##SBATCH --dependency=afterok:45514509

module load vasp/6.2.1-gpu

source ~/.bash_aliases; keep_log

srun -n 16 -c 32 --cpu-bind=cores -G 16 --gpu-bind=single:1 vasp_ncl
"""

job_script = job_script_perlmutter_4_nodes

emats_file = []

for i in range(100):
    emats_file.append( [[1,0,0],[0,1,0],[0,0,1]] )

emat_file = [ \
    [ 1, 0, 0 ], \
    [ 0, 1, 0 ], \
    [ 0, 0, 1 ] ] 

