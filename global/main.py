import os
import numpy as np
from common import *

e_ref = -1318.22495360

def find_global_directions():
    direction_emin = (4.92, 0)
    myjob = vasp_jobs_ncl(n_theta=1801, n_phi=3600)
    myjob.setup_local_ref_frames(z_direction=direction_emin, from_file=False)
    myjob.add_directions([[0.81,270]], local_ref_frame=True)
    myjob.print_dirs_and_configs(local_ref_frame=False)
    #print("Opposite direction: {:6.1f} {:6.1f}".format(*opposite_direction))

def do():

    direction_emin = (5.0, 350.6)

    myjob = vasp_jobs_ncl(n_theta=1801, n_phi=3600)
    myjob.setup_local_ref_frames(z_direction=direction_emin, from_file=False)
    #myjob.add_directions([[0,0]], local_ref_frame=False)
    #myjob.add_thetas(phi=180, theta_min=6, theta_max=0, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas(  phi=0, theta_min=0, theta_max=6, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas(phi=270, theta_min=6, theta_max=0, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas( phi=90, theta_min=0, theta_max=6, ntheta=4, local_ref_frame=True)
    myjob.add_thetas(phi=0, theta_min=0, theta_max=180, ntheta=7, local_ref_frame=True)
    myjob.add_thetas(phi=180, theta_min=0, theta_max=180, ntheta=7, local_ref_frame=True)
    myjob.add_thetas(phi=90, theta_min=0, theta_max=180, ntheta=7, local_ref_frame=True)
    myjob.add_thetas(phi=270, theta_min=0, theta_max=180, ntheta=7, local_ref_frame=True)
    myjob.add_phis(theta=90, phi_min=0, phi_max=330, nphi=12, local_ref_frame=True)

    myjob.e_ref = e_ref

    myjob.print_dirs_and_configs(local_ref_frame=True)
    #myjob.print_dirs_and_configs(local_ref_frame=False)
    #myjob.setup_jobs(submit=False)
    #myjob.setup_jobs(submit=True)
    #restart(myjob, test=True, max_energy=np.inf)
    #restart(myjob, test=False, max_energy=np.inf)
    #myjob.check_local_magmoms_all_configurations()
    #myjob.get_energies(max_energy=np.inf, de0=1.e-8)
    #get_all_energies(myjob, max_energy=np.inf, de0=1.e-7)

if __name__ == "__main__":

    #find_global_directions()

    do()

    pass

