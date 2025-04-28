import os
import numpy as np
from common import *

def find_global_directions():
    direction_emin = (71.4,  59.5)
    myjob = vasp_jobs_ncl(n_theta=1801, n_phi=3600)
    myjob.setup_local_ref_frames(z_direction=direction_emin, from_file=False)
    myjob.add_directions([[0.92, 180]], local_ref_frame=True)
    myjob.print_dirs_and_configs(local_ref_frame=False)
    #print("Opposite direction: {:6.1f} {:6.1f}".format(*opposite_direction))

def do(direction_emin = (0.0, 0.0), e_ref = 0.0, max_energy = np.inf, de0 = 1.e-8):

    myjob = vasp_jobs_ncl(n_theta=1801, n_phi=3600)
    #myjob.setup_local_ref_frames(z_direction=direction_emin, from_file=True)

    #myjob.add_directions([[0,0]], local_ref_frame=False)

    #myjob.add_thetas(  phi=0, theta_min=0, theta_max=30, ntheta=4, local_ref_frame=True)

    #myjob.add_thetas(phi=0, theta_min=0, theta_max=180, ntheta=7, local_ref_frame=True)
    #myjob.add_thetas(phi=180, theta_min=0, theta_max=180, ntheta=7, local_ref_frame=True)
    #myjob.add_thetas(phi=90, theta_min=0, theta_max=180, ntheta=7, local_ref_frame=True)
    #myjob.add_thetas(phi=270, theta_min=0, theta_max=180, ntheta=7, local_ref_frame=True)
    #myjob.add_phis(theta=90, phi_min=0, phi_max=330, nphi=12, local_ref_frame=True)

    #myjob.add_phis(theta=90, phi_min=315, phi_max=345, nphi=7, local_ref_frame=True)
    #myjob.add_thetas(phi=330, theta_min=75, theta_max=105, ntheta=7, local_ref_frame=True)

    myjob.setup_local_ref_frames(z_direction=direction_emin, from_file=False)
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

    #myjob.print_dirs_and_configs(local_ref_frame=True)
    #myjob.print_dirs_and_configs(local_ref_frame=False)
    #myjob.setup_jobs(submit=False)
    #myjob.setup_jobs(submit=True)
    #myjob.check_convergences(restart=False, de0=de0)
    #myjob.check_convergences(restart=True, de0=de0)
    #myjob.check_local_magmoms_all_configurations()
    #restart(myjob, test=True, max_energy=max_energy)
    #restart(myjob, test=False, max_energy=max_energy)
    myjob.get_energies(max_energy=max_energy, de0=de0)
    #get_all_energies(myjob, max_energy=max_energy, de0=de0)

if __name__ == "__main__":

    #find_global_directions()

    do(direction_emin = (98.282, 330.142), e_ref = -1834.36333292, max_energy = np.inf, de0 = 1.e-8)

    pass

