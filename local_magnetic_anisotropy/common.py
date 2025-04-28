import os
import sys
import subprocess
import math
import time
import numpy as np
import copy
from base import magmoms, sphere, change_frame_sph, get_opposite_direction_sph, get_angle, get_angle_sph, sph2cart, get_emat_local
from data import root_dir, incar, poscar, kpoints, job_script, emats_file

class vasp_job:

    def __init__(self, dir_name=""):
        # root_dir etc are imported from data.py
        self.root_dir = root_dir
        self.dir_name = dir_name
        self.incar = incar
        self.poscar = poscar
        self.kpoints = kpoints
        self.job_script = job_script
        self.e_ref = 0

    def set_dir_name(self, dir_name=""):
        self.dir_name = dir_name

    def create_dir(self):
        if not os.path.isdir(self.dir_name):
            subprocess.run(["mkdir", "-p", self.dir_name])

    def get_input_files(self, mstring="3*"):
        os.chdir(self.dir_name)
        subprocess.run(["pwd"])
        subprocess.run(["ln", "-sf", self.root_dir + "POTCAR", "."])
        with open("INCAR", "w") as f:
            f.write(incar.format(mstring = mstring, rwigs = self.get_rwigs()))
        with open("POSCAR", "w") as f:
            f.write(poscar)
        with open("KPOINTS", "w") as f:
            f.write(kpoints)
        with open("vasp.job", "w") as f:
            f.write(job_script)
        #subprocess.run(["ln", "-sf", "/global/common/software/nersc/pm-stable/sw/vasp/vdw_kernal/vdw_kernel.bindat", "."])
        os.chdir(self.root_dir)

    def get_rwigs(self):
        # Get the Wigner-Seitz radii of the atoms from POTCAR
        rwigs = ""
        try:
            out = subprocess.run(["grep", "RWIGS", "POTCAR"], capture_output=True, text=True)
            out = out.stdout.rstrip().split("\n")
            for i in range(len(out)):
                rwigs = "{:s} {:s}".format(rwigs, out[i].split()[5])
        except:
            print("Error in getting RWIGS from POTCAR. Stopping ...")
            exit()
        return rwigs

    def submit_job(self):
        os.chdir(self.dir_name)
        subprocess.run(["pwd"])
        subprocess.run(["sbatch", "vasp.job"])
        time.sleep(1)
        os.chdir(self.root_dir)

    def check_convergence(self, restart=False, de0=1.e-8):
        self.convergence = False
        os.chdir(self.dir_name)
        #subprocess.run(["pwd"])
        try:
            out = subprocess.run(["grep", ":", "output"], capture_output=True)
            out = out.stdout.decode("utf-8").split('\n')
            de = out[-2].split()[3]
            if abs(float(de)) < de0:
                self.convergence = True
            print("dir_name = {:>50s} , de = {:>15s} , convergence = {:d}.".format(self.dir_name, de, self.convergence))
        except:
            print("dir_name = {:>50s} , de = {:>15s} , convergence = {:d}.".format(self.dir_name, "NA", self.convergence))
            pass
        if (not self.convergence) and restart:
            subprocess.run(["sbatch", "vasp.job"])
        os.chdir(self.root_dir)

    def count_nstep(self):
        os.chdir(self.dir_name)
        #subprocess.run(["pwd"])
        self.nstep = 0
        try:
            out = subprocess.run(["grep", "DAV:", "output"], capture_output=True)
            steps = out.stdout.split(b'\n')
            nstep = len(steps)
            de = steps[-2].decode("utf-8").split()[3]
            if abs(float(de)) <= 1.e-8:
                self.nstep = nstep
        except:
            pass
        os.chdir(self.root_dir)

    def get_energy(self, de0=1.e-8, max_energy=0.01):
        self.energy = 0
        os.chdir(self.dir_name)
        #subprocess.run(["pwd"])
        try:
            out = subprocess.run(["grep", ":", "output"], capture_output=True)
            out = out.stdout.decode("utf-8").split('\n')
            de = out[-2].split()[3]
            print(de)
            if abs(float(de)) < de0:
                try:
                    out = subprocess.run(["grep", "sigma", "OUTCAR"], capture_output=True)
                    e = float( out.stdout.decode("utf-8").split()[-1] )
                    if (e - self.e_ref) < max_energy:
                        self.energy = e
                except:
                    pass
        except:
            pass
        os.chdir(self.root_dir)

    def get_dipole(self, crystal=False):

        self.dipole = [np.nan, np.nan, np.nan, np.nan]
        os.chdir(self.dir_name)
        #subprocess.run(["pwd"])
        if crystal:
            tag = "electronic dipole moment"
            m = 5
            n = 8
        else:
            tag = "dipolmoment"
            m = 1
            n = 4
        try:
            with open("OUTCAR", "r") as f:
                for line in f:
                    if tag in line:
                        tmp = line.strip().split()
                        self.dipole = [float(tmp[i]) for i in range(m, n)]
                        self.dipole.append(np.linalg.norm(self.dipole))
                        break
        except:
            pass
        os.chdir(self.root_dir)

class vasp_job_ncl(vasp_job):
    def __init__(self, n_site=3, magnetic_moments=[3, 3, 3], n_theta=1801, n_phi=3600):
        super().__init__()

        # self.mms and self.sph are in the global reference frame
        self.mms = magmoms(n_site=n_site, magmoms=magnetic_moments)
        self.sph = sphere(n_theta=n_theta, n_phi=n_phi)

        # example self.mms.directions_of_spin: [[0,0], [0,0], [0,0]]
        self.dir_name = self.get_dir_name(self.mms.directions_of_spin)

        # set self.z_directions, self.emats, and self.directions_of_spin_local
        self.setup_local_ref_frames()

    def get_dir_name(self, directions_of_spin):

        dir_name = ""
        for i in range(self.mms.n_site):
            self.sph.set_angle(directions_of_spin[i])
            dir_name = dir_name + str(self.sph.i_direction)
            if i < self.mms.n_site - 1:
                dir_name = dir_name + "_"
        return dir_name

    def set_dir_name(self, dir_name="0_0_0"):

        self.dir_name = dir_name

        i_directions = [int(dir_name.split("_")[i]) for i in range(self.mms.n_site)]
        directions_of_spin = []
        for i in range(self.mms.n_site):
            self.sph.set_i_direction(i_direction = i_directions[i])
            directions_of_spin.append( [self.sph.theta, self.sph.phi] )
        self.mms.set_directions_of_spin(directions_of_spin)

    def setup_local_ref_frames(self, sites=[], z_directions=[], from_file=False):
        # example input: (sites=[0], z_directions=[(0,0)], from_file=False)
        # Let the sphere of local reference frame be continous.
        # Let the sphere of global reference frame be discrete.

        # number of sites for which we specify local reference frames.
        n_site = len(sites)

        self.z_directions = [(0, 0) for i in range(self.mms.n_site)]
        self.emats = [np.eye(3) for i in range(self.mms.n_site)]

        for i in range(n_site):
            if from_file:
                self.emats[sites[i]] = emats_file[sites[i]]
            else:
                self.z_directions[sites[i]] = z_directions[i]
                self.emats[sites[i]] = get_emat_local(z_direction = z_directions[i])

        self.transform_directions_of_spin(global2local=True)

    def transform_directions_of_spin(self, global2local=True):
        # The local reference frames are defined by self.emats

        if global2local:
            self.directions_of_spin_local = []
            for i in range(self.mms.n_site):
                theta, phi = self.mms.directions_of_spin[i]
                dumb, theta, phi = change_frame_sph([1, theta, phi], np.eye(3), self.emats[i])
                self.directions_of_spin_local.append([theta, phi])
        else:
            directions_of_spin = []
            for i in range(self.mms.n_site):
                theta, phi = self.directions_of_spin_local[i]
                dumb, theta, phi = change_frame_sph([1, theta, phi], self.emats[i], np.eye(3))
                directions_of_spin.append([theta, phi])
            self.mms.set_directions_of_spin(directions_of_spin)

    def set_directions_of_spin(self, directions_of_spin=[[0,0],[0,0],[0,0]], local_ref_frame=False):
        # sets self.mms and self.directions_of_spin_local
        if local_ref_frame:
            self.directions_of_spin_local = directions_of_spin
            self.transform_directions_of_spin(global2local=False)
        else:
            self.mms.set_directions_of_spin(directions_of_spin)
            self.transform_directions_of_spin(global2local=True)

    def setup_one_job(self, directions_of_spin = [[0,0], [0,0], [0,0]], local_ref_frame=False, submit=False):

        self.set_directions_of_spin(directions_of_spin=directions_of_spin, local_ref_frame=local_ref_frame)
        self.dir_name = self.get_dir_name(self.mms.directions_of_spin)
        self.create_dir()
        if submit:
            self.submit_job()
        else:
            self.get_input_files(mstring = self.mms.mstring)

    def get_local_magmoms_ncl(n_site):
        f = open('OUTCAR', 'r')
        data = f.readlines()
        f.close()

        iline_mx = 0
        iline_my = 0
        iline_mz = 0

        for i, line in enumerate(data):
            if "magnetization (x)" in line:
                iline_mx = i
            if "magnetization (y)" in line:
                iline_my = i
            if "magnetization (z)" in line:
                iline_mz = i
                break

        magnetic_moments = []

        for i_site in range(n_site):
            magnetic_moments.append([])
            magnetic_moments[-1].append(float(data[iline_mx+4+i_site].strip().split()[-1]))
            magnetic_moments[-1].append(float(data[iline_my+4+i_site].strip().split()[-1]))
            magnetic_moments[-1].append(float(data[iline_mz+4+i_site].strip().split()[-1]))

        with open("local_magmoms", "w") as f:
            for i_site in range(n_site):
                f.write("{:6.2f} {:6.2f} {:6.2f}\n".format(*magnetic_moments[i_site]))

        return np.array( magnetic_moments )

    def check_local_magmoms(self):
        os.chdir(self.dir_name)

        #subprocess.run(["pwd"])

        try:
            mms = self.get_local_magmoms_ncl(self.mms.n_site)
            #print(mms)
            
            mms_init = []
            for i in range(self.mms.n_site):
                direction = self.mms.directions[i]
                mms_init.append( sph2cart([self.mms.magmoms[i], direction[0], direction[1]]) )
            mms_init = np.array(mms_init)
            #print(mms_init)
            
            angles = []
            for i in range(self.mms.n_site):
                angles.append( get_angle(mms[i], mms_init[i]) )
            print(("{:50s} " + "{:6.1f} "*self.mms.n_site).format(self.dir_name, *angles))
        except:
            print(("{:50s} " + "   Nan "*self.mms.n_site).format(self.dir_name))

        os.chdir(root_dir)

    def get_local_magmom(self, i_site):
        os.chdir(self.dir_name)

        try:
            f = open('OUTCAR', 'r')
            data = f.readlines()
            f.close()
            
            iline_mx = 0
            iline_my = 0
            iline_mz = 0
            
            for i, line in enumerate(data):
                if "magnetization (x)" in line:
                    iline_mx = i
                if "magnetization (y)" in line:
                    iline_my = i
                if "magnetization (z)" in line:
                    iline_mz = i
                    break
            
            lmm = []
            lmm.append(float(data[iline_mx+4+i_site].strip().split()[-1]))
            lmm.append(float(data[iline_my+4+i_site].strip().split()[-1]))
            lmm.append(float(data[iline_mz+4+i_site].strip().split()[-1]))
            lmm = np.array(lmm)
        except:
            lmm = np.array([np.nan, np.nan, np.nan])
        self.lmm = lmm

        os.chdir(root_dir)

    def get_local_orbmom(self, i_site):
        os.chdir(self.dir_name)

        try:
            f = open('OUTCAR', 'r')
            data = f.readlines()
            f.close()
            
            iline_mx = 0
            iline_my = 0
            iline_mz = 0
            
            for i, line in enumerate(data):
                if "orbital moment (x)" in line:
                    iline_mx = i
                if "orbital moment (y)" in line:
                    iline_my = i
                if "orbital moment (z)" in line:
                    iline_mz = i
                    break
            
            lom = []
            lom.append(float(data[iline_mx+4+i_site].strip().split()[-1]))
            lom.append(float(data[iline_my+4+i_site].strip().split()[-1]))
            lom.append(float(data[iline_mz+4+i_site].strip().split()[-1]))
            lom = np.array(lom)
        except:
            lom = np.array([np.nan, np.nan, np.nan])
        self.lom = lom

        os.chdir(root_dir)

    def get_local_orbmoms_ncl(n_site):
        os.chdir(self.dir_name)

        #subprocess.run(["pwd"])

        f = open('OUTCAR', 'r')
        data = f.readlines()
        f.close()

        iline_mx = 0
        iline_my = 0
        iline_mz = 0

        for i, line in enumerate(data):
            if "orbital moment (x)" in line:
                iline_mx = i
            if "orbital moment (y)" in line:
                iline_my = i
            if "orbital moment (z)" in line:
                iline_mz = i
                break

        orbital_moments = []

        for i_site in range(n_site):
            orbital_moments.append([])
            orbital_moments[-1].append(float(data[iline_mx+4+i_site].strip().split()[-1]))
            orbital_moments[-1].append(float(data[iline_my+4+i_site].strip().split()[-1]))
            orbital_moments[-1].append(float(data[iline_mz+4+i_site].strip().split()[-1]))

        with open("local_orbmoms.dat", "w") as f:
            for i_site in range(n_site):
                f.write("{:6.2f} {:6.2f} {:6.2f}\n".format(*orbital_moments[i_site]))

        print("Local orbital moments are written to local_orbmoms.dat")

        os.chdir(root_dir)

    def get_occmat():
        os.chdir(self.dir_name)

        subprocess.run(["pwd"])

        f = open('OUTCAR', 'r')
        data = f.readlines()
        f.close()

        iline_start = 0
        for i, line in enumerate(data):
            if 'Iteration' in line:
                iline_start = i

        data = data[iline_start:]

        natom = 0
        line_numbers = []
        ls = []
        for i, line in enumerate(data):
            if "atom =" in line and "l =" in line:
                ls.append(int(line.strip().split()[-1]))
                natom += 1
                line_numbers.append(i)

        skip = 4
        n_spin = 4 # For noncollinear DFT calculations

        entry = "# {:3d}\n".format(natom)
        i_atom = 0
        for i_atom in range(natom):
            l = ls[i_atom]
            tag = "# {:3d} {:3d} {:3d}\n".format(i_atom+1, l, n_spin)
            n_line = n_spin*(2*l+1 + 3)
            entry = entry + tag
            for line in data[line_numbers[i_atom]+skip:line_numbers[i_atom]+skip+n_line]:
                if line.strip() == "":
                    continue
                else:
                    entry += line
            entry += "# \n"

        entry = entry.replace('spin', '# spin')

        with open('occmat', 'w') as f:
            f.write(entry)

        print("Occupancy matrix is written to occmat")

        os.chdir(root_dir)

    def get_occmat_eigen():
        os.chdir(self.dir_name)

        subprocess.run(["pwd"])

        f = open('OUTCAR', 'r')
        data = f.readlines()
        f.close()

        iline_start = 0
        for i, line in enumerate(data):
            if 'Iteration' in line:
                iline_start = i

        data = data[iline_start:]

        natom = 0
        line_numbers = []
        ls = []
        for i, line in enumerate(data):
            if "atom =" in line and "l =" in line:
                ls.append(int(line.strip().split()[-1]))
                natom += 1
                line_numbers.append(i)

        skip = 38 # For noncollinear DFT calculations
        n_spin = 4 # For noncollinear DFT calculations
        i_atom = 0
        entry = ""
        for i_atom in range(natom):
            l = ls[i_atom]
            tag = "# {:3d} {:3d} {:3d}\n".format(i_atom+1, l, n_spin)
            n_line = 2*(2*l + 1)
            entry = entry + tag
            for line in data[line_numbers[i_atom]+skip:line_numbers[i_atom]+skip+n_line]:
                if line.strip() == "":
                    continue
                else:
                    line = line.replace("o =", "")
                    line = line.strip().split("v =")
                    line = line[1] + "   #   " + line[0] + "\n"
                    entry += line
            print( entry )

        with open('occmat.eigen', 'w') as f:
            f.write(entry)

        print("Eigenvalues and eigenvectors of occupancy matrix are written to occmat.eigen")

        os.chdir(root_dir)

class vasp_jobs_ncl(vasp_job_ncl):
    def __init__(self, n_site=3, magnetic_moments=[3, 3, 3], n_theta=1801, n_phi=3600):
        super().__init__(n_site=n_site, magnetic_moments=magnetic_moments, n_theta=n_theta, n_phi=n_phi)

        self.configurations = []
        self.configurations_local = []
        self.n_conf = len(self.configurations)
        self.dir_names = []
        self.energies = []
        self.dipoles = []

    def get_dir_names(self):
        dir_names = []
        for i in range(self.n_conf):
            directions_of_spin = self.configurations[i]
            dir_name = self.get_dir_name(directions_of_spin)
            dir_names.append(dir_name)
        self.dir_names = dir_names

    def add_configurations(self, configurations=[[[0,0],[0,0],[0,0]]], local_ref_frame=False):

        # All configurations share the same set of local reference frames
        # The default local reference frames are the same with the global one
        # The local reference frames can be set by self.setup_local_ref_frames

        n_conf = len( configurations )

        for i in range(n_conf):
            self.set_directions_of_spin(directions_of_spin=configurations[i], local_ref_frame=local_ref_frame)
            dir_name = self.get_dir_name(self.mms.directions_of_spin)
            if dir_name not in self.dir_names:
                self.n_conf = self.n_conf + 1
                self.configurations.append(self.mms.directions_of_spin)
                self.configurations_local.append(self.directions_of_spin_local)
                self.dir_names.append(dir_name)
                self.energies.append(0)
                self.dipoles.append([np.nan, np.nan, np.nan])

    def add_one_configuration_by_dir_name(self, dir_name="0_0_0"):

        self.set_dir_name(dir_name)
        self.add_configurations(configurations=[self.mms.directions_of_spin], local_ref_frame=False)

    def add_thetas_colinear_spin(self, phi=0, theta_min=0, theta_max=180, ntheta=2, local_ref_frame=False):
        thetas = np.linspace(theta_min, theta_max, ntheta, endpoint=True)
    
        for theta in thetas:
            directions_of_spin = [[theta, phi] for i in range(self.mms.n_site)] 
            self.set_directions_of_spin(directions_of_spin=directions_of_spin, local_ref_frame=local_ref_frame)
            dir_name = self.get_dir_name(self.mms.directions_of_spin)
            if dir_name not in self.dir_names:
                self.n_conf = self.n_conf + 1
                self.configurations.append(self.mms.directions_of_spin)
                self.configurations_local.append(self.directions_of_spin_local)
                self.dir_names.append(dir_name)
                self.energies.append(0)
                self.dipoles.append([np.nan, np.nan, np.nan])

    def add_phis_colinear_spin(self, theta=0, phi_min=0, phi_max=360, nphi=2, local_ref_frame=False):
        phis = np.linspace(phi_min, phi_max, nphi, endpoint=True)
    
        for phi in phis:
            directions_of_spin = [[theta, phi] for i in range(self.mms.n_site)] 
            self.set_directions_of_spin(directions_of_spin=directions_of_spin, local_ref_frame=local_ref_frame)
            dir_name = self.get_dir_name(self.mms.directions_of_spin)
            if dir_name not in self.dir_names:
                self.n_conf = self.n_conf + 1
                self.configurations.append(self.mms.directions_of_spin)
                self.configurations_local.append(self.directions_of_spin_local)
                self.dir_names.append(dir_name)
                self.energies.append(0)
                self.dipoles.append([np.nan, np.nan, np.nan])

    def set_direction_for_one_site(self, i_site=0, direction=[0,0], local_ref_frame=False):
        self.transform_configurations(global2local = True)
        if local_ref_frame:
            for i in range(self.n_conf):
                self.configurations_local[i][i_site] = direction
        else:
            for i in range(self.n_conf):
                self.configurations[i][i_site] = direction
        self.transform_configurations(global2local = (not local_ref_frame))
        self.get_dir_names()

    def flip_direction_for_one_site(self, i_site=0):
        for i in range(self.n_conf):
            theta, phi = self.configurations[i][i_site]
            theta, phi = get_opposite_direction_sph((theta, phi))
            self.configurations[i][i_site] = [theta, phi]
        self.transform_configurations(global2local=True)
        self.get_dir_names()

    def setup_jobs(self, submit=False):
        for i in range(self.n_conf):
            self.setup_one_job(directions_of_spin = self.configurations[i], local_ref_frame=False, submit=submit)

    def check_convergences(self, restart=False, de0=1.e-8):
        self.convergences = []
        for i in range(self.n_conf):
            self.dir_name = self.dir_names[i]
            self.set_directions_of_spin(directions_of_spin=self.configurations[i], local_ref_frame=False)
            self.check_convergence(restart=restart, de0=de0)
            self.convergences.append( self.convergence )

    def get_energies(self, max_energy=0.01, de0=1.e-8):
        self.energies = []

        ostring1 = ""
        ostring2 = ""
        for i in range(self.n_conf):
            print(self.dir_names[i])
            self.dir_name = self.dir_names[i]
            self.get_energy(de0=de0, max_energy=max_energy)
            self.energies.append(self.energy)
            ostring1 = ostring1 + "#"

            for j in range(self.mms.n_site):
                ostring1 = ostring1 + "{:6.1f}  {:6.1f}  ".format(self.configurations[i][j][0], self.configurations[i][j][1])
                ostring2 = ostring2 + "{:6.1f}  {:6.1f}  ".format(self.configurations_local[i][j][0], self.configurations_local[i][j][1])
            ostring1 = ostring1 + "{:15.9f}\n".format( self.energy )
            ostring2 = ostring2 + "{:15.9f}\n".format( self.energy )
    
        with open("energies.dat", "w") as f:
            f.write("# global reference frame\n")
            f.write(ostring1)
            f.write("\n\n")
            f.write("# local reference frame\n")
            f.write(ostring2)

    def get_dipoles(self, crystal=False):
        self.dipoles = []

        ostring1 = ""
        ostring2 = ""
        for i in range(self.n_conf):
            print(self.dir_names[i])
            self.dir_name = self.dir_names[i]
            self.get_dipole(crystal=crystal)
            self.dipoles.append(self.dipole)
            ostring1 = ostring1 + "#"

            for j in range(self.mms.n_site):
                ostring1 = ostring1 + "{:6.1f}  {:6.1f}  ".format(self.configurations[i][j][0], self.configurations[i][j][1])
                ostring2 = ostring2 + "{:6.1f}  {:6.1f}  ".format(self.configurations_local[i][j][0], self.configurations_local[i][j][1])
            ostring1 = ostring1 + "{:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format( *self.dipole)
            ostring2 = ostring2 + "{:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format( *self.dipole)
    
        with open("dipoles.dat", "w") as f:
            f.write("# global reference frame\n")
            f.write(ostring1)
            f.write("\n\n")
            f.write("# local reference frame\n")
            f.write(ostring2)

    def print_dirs_and_configs(self, local_ref_frame=False):
        # The global reference frame is emat0 = np.eye(3). 
        # The local reference frame is defined by self.emats.

        ostring = ""
        for i in range(self.mms.n_site):
            ostring = ostring + "{:6.1f} {:6.1f} "
        ostring = ostring + ": {:30s}"

        if local_ref_frame:
            configurations = self.configurations_local
        else:
            configurations = self.configurations

        for i in range(self.n_conf):
            directions_of_spin_flat = [item for direction in configurations[i] for item in direction]
            print(ostring.format(*directions_of_spin_flat, self.dir_names[i]))

    def transform_configurations(self, global2local=True):

        # Assume that self.emats (definition of local reference frames) are already specified.

        if global2local:
            configurations = self.configurations
        else:
            configurations = self.configurations_local

        for i in range(self.n_conf):
            self.set_directions_of_spin(directions_of_spin=configurations[i], local_ref_frame=(not global2local))
            self.configurations[i] = self.mms.directions_of_spin
            self.configurations_local[i] = self.directions_of_spin_local

    def check_local_magmoms_all_configurations(self):
        for i in range(self.n_conf):
            self.set_directions_of_spin(directions_of_spin=self.configurations[i], local_ref_frame=False)
            self.dir_name = self.dir_names[i]
            self.check_local_magmoms()

    def get_local_magmom_and_orbmom_all_configurations(self, i_site):
        self.lmms = []
        self.loms = []
        ostring1 = ""
        ostring2 = ""
        for i in range(self.n_conf):
            print(self.dir_names[i])
            self.set_directions_of_spin(directions_of_spin=self.configurations[i], local_ref_frame=False)
            self.dir_name = self.dir_names[i]
            self.get_local_magmom(i_site)
            self.get_local_orbmom(i_site)
            self.lmms.append(self.lmm)
            self.loms.append(self.lom)
            try:
                angle = get_angle(self.lmm, self.lom, rad=False)
            except:
                angle = np.nan

            ostring1 = ostring1 + "#"

            for j in range(self.mms.n_site):
                ostring1 = ostring1 + "{:6.1f}  {:6.1f}  ".format(self.configurations[i][j][0], self.configurations[i][j][1])
                ostring2 = ostring2 + "{:6.1f}  {:6.1f}  ".format(self.configurations_local[i][j][0], self.configurations_local[i][j][1])
            ostring1 = ostring1 + "{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format( *self.lmm, *self.lom, angle)
            ostring2 = ostring2 + "{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format( *self.lmm, *self.lom, angle)
    
        with open("lmms_and_loms_site" + str(i_site) + ".dat", "w") as f:
            f.write("# global reference frame\n")
            f.write(ostring1)
            f.write("\n\n")
            f.write("# local reference frame\n")
            f.write(ostring2)

    def get_occmat_eigenvectors_all_configurations(self):
        for i in range(self.n_conf):
            self.dir_name = self.dir_names[i]
            self.get_occmat_eigenvectors()

def restart(myjob, test=True, max_angle=180, de0=1.e-8, from_neighbor=True, file_name="CHGCAR", copy_it=True):
    # file_name = "WAVECAR" or "CHGCAR"

    myjob.check_convergences(restart=False, de0=de0)

    jobs_already_done = [i for i, x in enumerate(myjob.convergences) if x == True]
    jobs_not_done = [i for i, x in enumerate(myjob.convergences) if x == False]

    if from_neighbor:
    
        ostring = "Restart ("
        unit_of_ostring = "({:5.1f},{:5.1f}),"
        for i in range(myjob.mms.n_site-1):
            ostring = ostring + unit_of_ostring
        ostring = ostring + "({:5.1f},{:5.1f})) / {:30s} from ("
        for i in range(myjob.mms.n_site-1):
            ostring = ostring + unit_of_ostring
        ostring = ostring + "({:5.1f},{:5.1f})) / {:30s}"
    
        for i in jobs_not_done:
            max_angles = []
            for j in jobs_already_done:
                angles = [float("nan") for i in range(myjob.mms.n_site)]
                for k in range(myjob.mms.n_site):
                    v1 = [1, myjob.configurations[i][k][0], myjob.configurations[i][k][1]]
                    v2 = [1, myjob.configurations[j][k][0], myjob.configurations[j][k][1]]
                    angles[k] = get_angle_sph(v1, v2)
                tmp = max(angles)
                max_angles.append((j, tmp))
            dtype = [("index", int), ("angle", float)]
            max_angles = np.array(max_angles, dtype=dtype)
            max_angles = np.sort(max_angles, order="angle")
            j = max_angles[0][0]
            angle = max_angles[0][1]
            #if test:
                #print(i, angle)
            if angle <= max_angle:
                directions_flat_i = [item for direction in myjob.configurations_local[i] for item in direction]
                directions_flat_j = [item for direction in myjob.configurations_local[j] for item in direction]
                orig = root_dir + myjob.dir_names[j] + "/" + file_name
                dest = root_dir + myjob.dir_names[i] + "/" + file_name
                print(ostring.format(*directions_flat_i, myjob.dir_names[i], *directions_flat_j, myjob.dir_names[j]))
                if not test:
                    if copy_it:
                        subprocess.run(["cp", orig, dest])
                    else:
                        subprocess.run(["ln", "-sf", orig, dest])
                    myjob.set_dir_name(dir_name = myjob.dir_names[i])
                    myjob.submit_job()
    else:
        ostring = "Restart ("
        unit_of_ostring = "({:5.1f},{:5.1f}),"
        for i in range(myjob.mms.n_site-1):
            ostring = ostring + unit_of_ostring
        ostring = ostring + "({:5.1f},{:5.1f})) / {:30s}"

        for i in jobs_not_done:
            directions_flat_i = [item for direction in myjob.configurations_local[i] for item in direction]
            print(ostring.format(*directions_flat_i, myjob.dir_names[i]))
            if not test:
                myjob.set_dir_name(dir_name = myjob.dir_names[i])
                myjob.submit_job()

def get_all_energies(myjob, max_energy=0.1, de0=1.e-8):
    # Assume all directories beginning with 0-9 are DFT directories.

    dir_list = os.listdir(".")
    for i in range(len(dir_list)):
        if os.path.isdir(dir_list[i]) and dir_list[i].split("_")[0].isnumeric():
            myjob.add_one_configuration_by_dir_name(dir_list[i])
    myjob.print_dirs_and_configs(local_ref_frame=True)
    myjob.get_energies(max_energy, de0)

if __name__ == "__main__":

    pass

