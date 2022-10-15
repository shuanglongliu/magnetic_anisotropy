import os
import sys
import subprocess
import math
import time
import numpy as np
import copy
from ase.io import read, write
from base import sph2cart, cart2sph, change_frame_sph, get_emat_local, magmoms, sphere
from data import root_dir, incar, poscar, kpoints, job_script, emats_file, emat_file

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

    def get_input_files(self, sstring="0 0 1"):
        os.chdir(self.dir_name)
        subprocess.run(["pwd"])
        with open("INCAR", "w") as f:
            f.write(incar.format(sstring = sstring))
        with open("POSCAR", "w") as f:
            f.write(poscar)
        with open("KPOINTS", "w") as f:
            f.write(kpoints)
        with open("vasp.job", "w") as f:
            f.write(job_script)
        subprocess.run(["ln", "-sf", self.root_dir + "POTCAR", "."])
        subprocess.run(["ln", "-sf", self.root_dir + "0/WAVECAR", "."])
        os.chdir(self.root_dir)

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
            print("dir_name = {:>15s} , de = {:>15s} , convergence = {:d}.".format(self.dir_name, de, self.convergence))
            if (not self.convergence) and restart:
                subprocess.run(["sbatch", "vasp.job"])
        except:
            pass
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

    def get_energy(self, max_energy=0.01, de0=1.e-8):
        self.energy = 0
        os.chdir(self.dir_name)
        #subprocess.run(["pwd"])
        try:
            out = subprocess.run(["grep", ":", "output"], capture_output=True)
            out = out.stdout.decode("utf-8").split('\n')
            de = out[-2].split()[3]
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

class vasp_job_ncl(vasp_job):
    def __init__(self, saxis=[0, 0], n_theta=1801, n_phi=3600):
        super().__init__()

        # self.sph are in the global reference frame
        self.sph = sphere(n_theta=n_theta, n_phi=n_phi)

        self.saxis = saxis

        self.set_sstring(saxis)

        self.dir_name = self.get_dir_name(saxis)

        self.setup_local_ref_frames()

    def set_sstring(self, saxis):

        saxis_cart = sph2cart([1, saxis[0], saxis[1]])
        self.sstring = "{:10.6f} {:10.6f} {:10.6f}".format(*saxis_cart)

    def get_dir_name(self, saxis):

        self.sph.set_angle(saxis)
        dir_name = str(self.sph.i_direction)
        return dir_name

    def set_dir_name(self, dir_name="0"):

        self.dir_name = dir_name

        i_direction = int(dir_name)
        self.sph.set_i_direction(i_direction = i_direction)
        self.set_saxis( [self.sph.theta, self.sph.phi] )

    def setup_one_job(self, saxis = [0,0], submit=False):

        self.set_sstring(saxis)

        self.dir_name = self.get_dir_name(saxis)
        self.create_dir()
        if submit:
            self.submit_job()
        else:
            self.get_input_files(sstring = self.sstring)

    def setup_local_ref_frames(self, z_direction=[0, 0], from_file=False):
        # Let the sphere of local reference frame be continous.
        # Let the sphere of global reference frame be discrete.

        if from_file:
            self.emat = emat_file
            r_, theta, phi = cart2sph(emat_file[2])
            self.z_direction = [theta, phi]
        else:
            self.z_direction = z_direction
            self.emat = get_emat_local(z_direction = z_direction)

    def transform_saxis(self, global2local=True):
        # The local reference frames are defined by self.emat

        if global2local:
            self.saxis_local = []
            theta, phi = self.saxis
            dumb, theta, phi = change_frame_sph([1, theta, phi], np.eye(3), self.emat)
            self.saxis_local = [theta, phi]
        else:
            saxis = []
            theta, phi = self.saxis_local
            dumb, theta, phi = change_frame_sph([1, theta, phi], self.emat, np.eye(3))
            self.saxis = [theta, phi]

    def set_saxis(self, saxis=[0,0], local_ref_frame=False):
        if local_ref_frame:
            self.saxis_local = saxis
            self.transform_saxis(global2local=False)
        else:
            self.saxis = saxis
            self.transform_saxis(global2local=True)

        saxis_cart = sph2cart([1, self.saxis[0], self.saxis[1]])
        self.sstring = "{:10.6f} {:10.6f} {:10.6f}".format(*saxis_cart)

class vasp_jobs_ncl(vasp_job_ncl):
    def __init__(self, n_theta=1801, n_phi=3600):
        super().__init__(n_theta=n_theta, n_phi=n_phi)

        self.directions = []
        self.directions_local = []
        self.n_direction = len(self.directions)
        self.dir_names = []
        self.energies = []

    def get_dir_names(self):
        dir_names = []
        for i in range(self.n_direction):
            saxis = self.directions[i]
            dir_name = self.get_dir_name(saxis)
            dir_names.append(dir_name)
        self.dir_names = dir_names

    def add_directions(self, directions=[[0,0]], local_ref_frame=False):

        n_direction = len( directions )

        for i in range(n_direction):
            self.set_saxis(saxis=directions[i], local_ref_frame=local_ref_frame)
            dir_name = self.get_dir_name(self.saxis)
            if dir_name not in self.dir_names:
                self.directions.append(self.saxis)
                self.directions_local.append(self.saxis_local)
                self.n_direction = self.n_direction + 1
                self.dir_names.append(dir_name)
                self.energies.append(0)

    def add_thetas(self, phi=0, theta_min=0, theta_max=180, ntheta=2, local_ref_frame=False):
        thetas = np.linspace(theta_min, theta_max, ntheta, endpoint=True)
    
        for theta in thetas:
            self.set_saxis(saxis=[theta, phi], local_ref_frame=local_ref_frame)
            dir_name = self.get_dir_name(self.saxis)
            if dir_name not in self.dir_names:
                self.directions.append(self.saxis)
                self.directions_local.append(self.saxis_local)
                self.n_direction = self.n_direction + 1
                self.dir_names.append(dir_name)
                self.energies.append(0)

    def add_phis(self, theta=0, phi_min=0, phi_max=360, nphi=2, local_ref_frame=False):
        phis = np.linspace(phi_min, phi_max, nphi, endpoint=True)
    
        for phi in phis:
            self.set_saxis(saxis=[theta, phi], local_ref_frame=local_ref_frame)
            dir_name = self.get_dir_name(self.saxis)
            if dir_name not in self.dir_names:
                self.directions.append(self.saxis)
                self.directions_local.append(self.saxis_local)
                self.n_direction = self.n_direction + 1
                self.dir_names.append(dir_name)
                self.energies.append(0)

    def setup_jobs(self, submit=False):
        for i in range(self.n_direction):
            self.setup_one_job(saxis = self.directions[i], submit=submit)

    def check_convergences(self):
        self.convergences = []
        for i in range(self.n_direction):
            self.dir_name = self.dir_names[i]
            self.set_saxis(self.directions[i], local_ref_frame=False)
            self.check_convergence()
            self.convergences.append( self.convergence )

    def get_energies(self, max_energy=0.01, de0=1.e-6):
        self.energies = []

        ostring1 = ""
        ostring2 = ""
        for i in range(self.n_direction):
            print(self.dir_names[i])
            self.dir_name = self.dir_names[i]
            self.get_energy(max_energy=max_energy, de0=de0)
            self.energies.append(self.energy)
            ostring1 = ostring1 + "#{:6.1f}  {:6.1f}  {:15.9f}\n".format(self.directions[i][0], self.directions[i][1], self.energy)
            ostring2 = ostring2 + "{:6.1f}  {:6.1f}  {:15.9f}\n".format(self.directions_local[i][0], self.directions_local[i][1], self.energy)
    
        with open("energies.dat", "w") as f:
            f.write("# lab reference frame\n")
            f.write(ostring1)
            f.write("\n\n# local reference frame\n")
            f.write(ostring2)

    def print_dirs_and_configs(self, local_ref_frame=False):
        # The global reference frame is emat0 = np.eye(3). 
        # The local reference frame is defined by self.emat.

        ostring = "{:6.1f} {:6.1f} : {:30s}"

        if local_ref_frame:
            directions = self.directions_local
        else:
            directions = self.directions

        for i in range(self.n_direction):
            print(ostring.format(directions[i][0], directions[i][1],  self.dir_names[i]))

    def transform_directions(self, global2local=True):

        # Assume that self.emat (definition of local reference frame) is already specified.

        if global2local:
            directions = self.directions
        else:
            directions = self.directions_local

        for i in range(self.n_direction):
            self.set_saxis(saxis=directions[i], local_ref_frame=(not global2local))
            self.directions[i] = self.saxis
            self.directions_local[i] = self.saxis_local

def restart(myjob, test=True, max_energy=0.01):

    myjob.get_energies(max_energy=max_energy)

    jobs_status = [False for i in range(myjob.n_direction)] 
    for i in range(myjob.n_direction):
        e = myjob.energies[i]
        if e - myjob.e_ref < max_energy:
            jobs_status[i] = True

    if test:
        for i in range(myjob.n_direction):
            print("{:3d}   {:30s}   {:b}".format(i, myjob.dir_names[i], jobs_status[i]))

    jobs_already_done = [i for i, x in enumerate(jobs_status) if x == True]
    jobs_not_done = [i for i, x in enumerate(jobs_status) if x == False]
    
    ostring = "Restart ({:5.1f}, {:5.1f}) / {:30s} from ({:5.1f},{:5.1f}) / {:30s}"

    for i in jobs_not_done:
        orig = root_dir + "0/WAVECAR"
        dest = root_dir + myjob.dir_names[i] + "/"
        if test:
            print(ostring.format(myjob.directions[i][0], myjob.directions[i][1], myjob.dir_names[i], 0, 0, "0"))
            #print("ln -s", orig, dest)
        else:
            print(ostring.format(myjob.directions[i][0], myjob.directions[i][1], myjob.dir_names[i], 0, 0, "0"))
            #subprocess.run(["ln", "-s", orig, dest])
            myjob.set_dir_name(dir_name = myjob.dir_names[i])
            myjob.submit_job()

if __name__ == "__main__":

    myjob = vasp_job_ncl(saxis = [90, 0])

    pass

