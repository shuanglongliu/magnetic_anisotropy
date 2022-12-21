import os
import numpy as np
from common import *

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

def do(direction_emin = (0.0, 0.0), e_ref = 0.0, max_energy = np.inf, de0 = 1.e-8):

    myjob = vasp_jobs_ncl(n_site=3, magnetic_moments=[2, 2, 3], n_theta=1801, n_phi=3600)
    myjob.setup_local_ref_frames(sites=[0, 1, 2], from_file=True)

    ### Get initial CHGCAR or WAVECAR.

    #myjob.add_configurations([[[  0,0], [  0,0], [  0,0]]], local_ref_frame=False)
    #myjob.add_configurations([[[180,0], [180,0], [180,0]]], local_ref_frame=False)

    ### Search for the local z-axis.

    #myjob.setup_local_ref_frames(sites=[0, 1, 2], z_directions=[direction_emin, direction_emin, direction_emin], from_file=False)

    #myjob.add_thetas_colinear_spin(phi=180, theta_min=6, theta_max=0, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(  phi=0, theta_min=0, theta_max=6, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(phi=270, theta_min=6, theta_max=0, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin( phi=90, theta_min=0, theta_max=6, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(phi=180, theta_min=174, theta_max=180, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(  phi=0, theta_min=180, theta_max=174, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin(phi=270, theta_min=174, theta_max=180, ntheta=4, local_ref_frame=True)
    #myjob.add_thetas_colinear_spin( phi=90, theta_min=180, theta_max=174, ntheta=4, local_ref_frame=True)

    ### Sample three great circles.

    myjob.add_thetas_colinear_spin(  phi=0, theta_min= 0, theta_max=180, ntheta= 7, local_ref_frame=True)
    myjob.add_thetas_colinear_spin(phi=180, theta_min= 0, theta_max=180, ntheta= 7, local_ref_frame=True)
    myjob.add_thetas_colinear_spin( phi=90, theta_min= 0, theta_max=180, ntheta= 7, local_ref_frame=True)
    myjob.add_thetas_colinear_spin(phi=270, theta_min= 0, theta_max=180, ntheta= 7, local_ref_frame=True)
    myjob.add_phis_colinear_spin( theta=90,   phi_min= 0,   phi_max=330,   nphi=12, local_ref_frame=True)

    ### Exchange interaction.

    #myjob.add_configurations([[[90, 0], [180, 0], [180, 0], [ 90, 0]]], local_ref_frame=True)
    #myjob.add_configurations([[[90,90], [180, 0], [180, 0], [ 90, 0]]], local_ref_frame=True)
    #myjob.add_configurations([[[ 0, 0], [180, 0], [180, 0], [ 90, 0]]], local_ref_frame=True)
    #myjob.add_configurations([[[90, 0], [180, 0], [180, 0], [ 90,90]]], local_ref_frame=True)
    #myjob.add_configurations([[[90,90], [180, 0], [180, 0], [ 90,90]]], local_ref_frame=True)
    #myjob.add_configurations([[[ 0, 0], [180, 0], [180, 0], [ 90,90]]], local_ref_frame=True)
    #myjob.add_configurations([[[90, 0], [180, 0], [180, 0], [180, 0]]], local_ref_frame=True)
    #myjob.add_configurations([[[90,90], [180, 0], [180, 0], [180, 0]]], local_ref_frame=True)
    #myjob.add_configurations([[[ 0, 0], [180, 0], [180, 0], [180, 0]]], local_ref_frame=True)
    #myjob.flip_direction_for_one_site(i_site=0)                     
                                                                    
    #myjob.add_configurations([[[90, 0], [180, 0], [180, 0], [ 90, 0]]], local_ref_frame=True)
    #myjob.add_configurations([[[90,90], [180, 0], [180, 0], [ 90, 0]]], local_ref_frame=True)
    #myjob.add_configurations([[[ 0, 0], [180, 0], [180, 0], [ 90, 0]]], local_ref_frame=True)
    #myjob.add_configurations([[[90, 0], [180, 0], [180, 0], [ 90,90]]], local_ref_frame=True)
    #myjob.add_configurations([[[90,90], [180, 0], [180, 0], [ 90,90]]], local_ref_frame=True)
    #myjob.add_configurations([[[ 0, 0], [180, 0], [180, 0], [ 90,90]]], local_ref_frame=True)
    #myjob.add_configurations([[[90, 0], [180, 0], [180, 0], [180, 0]]], local_ref_frame=True)
    #myjob.add_configurations([[[90,90], [180, 0], [180, 0], [180, 0]]], local_ref_frame=True)
    #myjob.add_configurations([[[ 0, 0], [180, 0], [180, 0], [180, 0]]], local_ref_frame=True)

    ### Fix or flip one spin for all configurations.

    #myjob.set_direction_for_one_site(i_site=0, direction=[180, 0], local_ref_frame=True)
    #myjob.set_direction_for_one_site(i_site=1, direction=[180, 0], local_ref_frame=True)
    #myjob.set_direction_for_one_site(i_site=2, direction=[180, 0], local_ref_frame=True)
    #myjob.flip_direction_for_one_site(i_site=2)

    ### If more than one local reference frames are used, it is necessary to update the local coordinates with respect to the last local reference frame.

    #myjob.transform_configurations(global2local = True)

    ### Set reference energy for checking if the total energy is reasonable. e_ref is used by get_energies and restart.

    myjob.e_ref = e_ref

    ### Let's go! 

    #myjob.print_dirs_and_configs(local_ref_frame=True)
    #myjob.print_dirs_and_configs(local_ref_frame=False)
    #myjob.setup_jobs(submit=False)
    #myjob.setup_jobs(submit=True)
    #myjob.check_convergences(restart=False, de0=de0)
    #myjob.check_convergences(restart=True, de0=de0)
    #myjob.check_local_magmoms_all_configurations()
    #myjob.get_occmat_eigenvectors_all_configurations()
    #restart(myjob, test=True, max_angle=45, max_energy=max_energy, from_neighbor=True, file_name="WAVECAR", copy_it=False)
    #restart(myjob, test=False, max_angle=45, max_energy=max_energy, from_neighbor=True, file_name="WAVECAR", copy_it=False)
    myjob.get_energies(max_energy=max_energy, de0=de0)
    #get_all_energies(myjob, max_energy=max_energy, de0=de0)

if __name__ == "__main__":

    #check_and_resubmit()

    #find_global_directions()

    do(direction_emin = (0.0, 0.0), e_ref = 0.0, max_energy = np.inf, de0 = 1.e-8)

    pass

