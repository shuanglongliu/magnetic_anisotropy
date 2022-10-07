root_dir = "/global/cfs/cdirs/m3346/shlufl/IETS/nco_basecell_beta_strain/mae/"

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

poscar = """Ni Co  O                                
    1.000000000000000     
     4.0445000000000002   -4.0445000000000002    0.0000000000000000
     4.0445000000000002    4.0445000000000002    0.0000000000000000
     0.0000000000000000    0.0000000000000000    8.1140000000000008
   Ni   Co   O 
     4     8    16
Direct
  0.2499994644141950  0.7499997086283940  0.2499988093728049
  0.7499987767811902  0.7499993972013286  0.2499994421086882
  0.7499988271954976  0.2499998627078526  0.7499998713739515
  0.2500025323587778  0.2499996938790190  0.7499994063103941
  0.9999977600593937  0.2500000558934588  0.1303674027635253
  0.5000019233469359  0.2500005482595640  0.3696387799268805
  0.5000024154929150  0.7500009425258582  0.6303619779112140
  0.9999975357320992  0.7500005731104977  0.8696375686067768
  0.0000018024598631  0.0000001787876371  0.4999998996470865
  0.0000019275282597  0.4999999967247888  0.5000006927687366
  0.4999980708497844  0.4999997664342573  0.9999995634558303
  0.4999982595545092  0.0000007537292603  0.9999998344024945
  0.2742804645773376  0.7500000841326298  0.0062246557675394
  0.7257167486464624  0.7499998848056464  0.0062258775100261
  0.4999978648907515  0.5156172487063841  0.2355221687682700
  0.4999995194468241  0.9843827929590390  0.2355222807207795
  0.9999991518378550  0.5156177299272997  0.2644761079124081
  0.0000001846151392  0.9843814618180389  0.2644751654455817
  0.7742879988042688  0.7499999899869110  0.4937770684404370
  0.2257169591519244  0.7500000951663566  0.4937765715735338
  0.7742851743660566  0.2499997695133374  0.5062245743929168
  0.2257183276672450  0.2499998652392676  0.5062244602962949
  0.0000000340953719  0.0156149177504119  0.7355215362918486
  0.0000018081668003  0.4843833419755015  0.7355219927486587
  0.4999995058645084  0.0156148199654993  0.7644787321401054
  0.5000010588142558  0.4843859348775581  0.7644781014493347
  0.2742795658727175  0.2500001417142030  0.9937736797192045
  0.7257163374090538  0.2500004435799923  0.9937737781746634
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

module purge; module load cuda/11.1.0 nvhpc/20.11 openmpi/4.0.5 qd/2.3.22 fftw/3.3.8 vasp/6.2.0

source ~/.bash_aliases; keep_log

srun --mpi=pmix vasp_ncl
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

