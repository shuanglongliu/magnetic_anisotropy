import os
import numpy as np
from common import *

e_ref = -2277.46276603

def check_and_resubmit():

    dirnames = [ "4319101_0_0", "3238801_0_0" ]
    n = len(dirnames)
    myjob = vasp_job()
    for i in range(n):
        myjob.set_dir_name(dirnames[i])
        #print(myjob.dir_name)
        myjob.check_convergence(restart=False)
        #myjob.check_convergence(restart=True)

def find_global_directions():
    direction_emin = (42.9,  102.5)
    myjob = vasp_jobs_ncl(n_site=1, magnetic_moments=[3], n_theta=1801, n_phi=3600)
    myjob.setup_local_ref_frames(sites=[0], z_directions=[direction_emin], from_file=False)
    myjob.add_configurations([[[2,90]]], local_ref_frame=True)
    myjob.print_dirs_and_configs(local_ref_frame=False)

def do():

    myjob = vasp_jobs_ncl(n_site=3, magnetic_moments=[3, 4, 4], n_theta=1801, n_phi=3600)
    myjob.setup_local_ref_frames(sites=[0, 1, 2], from_file=True)

    myjob.add_configurations([[[  0,0], [  0,0], [  0,0]]], local_ref_frame=True)
    myjob.add_configurations([[[180,0], [180,0], [180,0]]], local_ref_frame=False)

    #direction_emin = (42.9, 102.4)
    #myjob.setup_local_ref_frames(sites=[0, 1, 2], z_directions=[direction_emin, direction_emin], from_file=False)
    #myjob.add_thetas_colinear_spin(phi=180, theta_min=6, theta_max=0, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(  phi=0, theta_min=0, theta_max=6, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(phi=270, theta_min=6, theta_max=0, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin( phi=90, theta_min=0, theta_max=6, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(phi=180, theta_min=174, theta_max=180, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(  phi=0, theta_min=180, theta_max=174, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(phi=270, theta_min=174, theta_max=180, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin( phi=90, theta_min=180, theta_max=174, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(  phi=0, theta_min= 0, theta_max=180, ntheta= 7, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(phi=180, theta_min= 0, theta_max=180, ntheta= 7, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin( phi=90, theta_min= 0, theta_max=180, ntheta= 7, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(phi=270, theta_min= 0, theta_max=180, ntheta= 7, local_ref_frame=True)
    #myjob.add_phis_colinear_spin( theta=90,   phi_min= 0,   phi_max=330,   nphi=12, local_ref_frame=True)

    #myjob.set_direction_for_one_site(i_site=1, direction=[0, 0], local_ref_frame=True)
    #myjob.set_direction_for_one_site(i_site=2, direction=[180, 0], local_ref_frame=True)
    #myjob.flip_direction_for_one_site(i_site=2)

    #myjob.transform_configurations(global2local = True)

    myjob.e_ref = e_ref

    #myjob.print_dirs_and_configs(local_ref_frame=True)
    #myjob.print_dirs_and_configs(local_ref_frame=False)
    #myjob.setup_jobs(submit=False)
    #myjob.setup_jobs(submit=True)
    #restart(myjob, test=True, max_angle=45, max_energy=0.1, from_neighbor=True, file_name="CHGCAR", copy_it=False)
    #restart(myjob, test=False, max_angle=45, max_energy=0.1, from_neighbor=True, file_name="CHGCAR", copy_it=False)
    #myjob.check_local_magmoms_all_configurations()
    #myjob.get_occmat_eigenvectors_all_configurations()
    myjob.get_energies(max_energy=0.1, de0=1.e-5)

if __name__ == "__main__":

    #check_and_resubmit()

    #find_global_directions()

    do()

    pass

