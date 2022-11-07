import os
import sys
import numpy as np
from scipy.optimize import curve_fit
from scipy.spatial.transform import Rotation as R
import matplotlib as mpl
import matplotlib.pyplot as plt
from spin import Spin, dim, Sx, Sy, Sz, ID, ss, Sx2, Sy2, Sz2
from pylib import change_frame_sph_rad, sph2cart_rad, cart2sph_deg

def rotate_frame(emat=np.eye(3), dtheta=0., dphi=0.):

    """
    Rotate the local reference frame
      rotate the z axis around the y axis by detheta, and then
      rotate the x and y axes around the new z axis by dphi
    Unit for angle: rad
    """

    # =================================================
    # rotate the z axis around the y axis by detheta
    # =================================================

    r = R.from_rotvec( dtheta * emat[1] )
    rotmat = r.as_matrix()
    ex = np.matmul(rotmat, emat[0])
    ey = emat[1]
    ez = np.matmul(rotmat, emat[2])

    # =======================================================
    # rotate the x and y axes around the new z axis by dphi
    # =======================================================

    r = R.from_rotvec( dphi * ez )
    rotmat = r.as_matrix()
    ex = np.matmul(rotmat, ex)
    ey = np.matmul(rotmat, ey)

    emat_prime = np.array( [ex, ey, ez] )

    return emat_prime


def change_frame(angles, emat=np.eye(3), dtheta=0., dphi=0.):

    """
    Unit for angle: rad
    """

    # ====================
    # Find new frame
    # ====================

    emat_prime = rotate_frame(emat, dtheta, dphi)

    # ====================
    # Change frame
    # ====================

    thetas, phis = angles

    thetas_prime = []
    phis_prime = []
    n_direction = len(thetas)
    for i in range(n_direction):
        theta = thetas[i]
        phi = phis[i]
        dumb, theta_prime, phi_prime = change_frame_sph_rad([1, theta, phi], emat, emat_prime)
        thetas_prime.append(theta_prime)
        phis_prime.append(phi_prime)
    thetas_prime = np.array(thetas_prime)
    phis_prime = np.array(phis_prime)

    return (thetas_prime, phis_prime)

def e_ani_q(angles, D, E, E0, dtheta, dphi):

    """
    Magnetic anisotropy energy by treating spins as quantum operators.
    Spins are treated as operators.
    Units for angle: rad
    """

    # ====================
    # Change frame
    # ====================

    thetas, phis = change_frame(angles, emat_dft, dtheta, dphi)


    # ==================================
    # Evaluate energy
    # ==================================

    energies = []
    n_direction = len(thetas)
    for i in range(n_direction):
        ex, ey, ez = sph2cart_rad([1, thetas[i], phis[i]])
    
        op = ex*Sx + ey*Sy + ez*Sz
    
        w, v = np.linalg.eig(op)
        w = np.real(w)
        #print((" {:12.6f}"*dim).format(*w))
        #print(v)
    
        indices = np.argsort(w)
        i = dim-1
        state = v[:, indices[i]]
        #print(w[indices[i]])
        #print(state); exit()

        op = 3*Sz2 + ss*ID
        energy = D/3*np.dot(np.conjugate(state), np.matmul(op, state))

        op = Sx2 - Sy2
        energy = energy + E*np.dot(np.conjugate(state), np.matmul(op, state))

        energy = energy + E0

        energies.append(energy)
    energies = np.real(energies)
    #print(D, energies); exit()

    return energies

def e_ani_c(directions, D, E, E0, dtheta, dphi):

    """
    Magnetic anisotropy energy by treating spins as classical vectors.
    Spins are treated as classical vectors.
    Units for angles: rad
    """

    global Spin, emat_dft

    thetas, phis = change_frame(angles, emat_dft, dtheta, dphi)

    return Spin**2 * ( D * np.cos(thetas)**2 + E * np.sin(thetas)**2 * np.cos(2*phis) ) + E0

def get_local_eigen_frame(dtheta, dphi):

    """
    Units for angle: rad
    """

    global emat_dft

    emat_loc = rotate_frame(emat_dft, dtheta, dphi)

    print("\nLocal eigen frame in spherical coordinates: ")
    dumb, theta, phi = cart2sph_deg(emat_loc[0])
    print("Local e_x: (theta, phi) = ({:6.1f}, {:6.1f}) deg.".format(theta, phi))
    dumb, theta, phi = cart2sph_deg(emat_loc[1])
    print("Local e_y: (theta, phi) = ({:6.1f}, {:6.1f}) deg.".format(theta, phi))
    dumb, theta, phi = cart2sph_deg(emat_loc[2])
    print("Local e_z: (theta, phi) = ({:6.1f}, {:6.1f}) deg.".format(theta, phi))

    print("\nLocal eigen frame in cartesian coordinates: ")
    tmp = np.concatenate(emat_loc).flat
    print((3*"{:12.6f} {:12.6f} {:12.6f}\n").format(*tmp))

    return emat_loc

def fit(base_name = "energies.dat", columns=[1, 2, 7], factor=1.0, title="OCHNH2, Co2", p0=[0, 0, 0, 0, 0], bounds=([0, 0, 0, 0, 0], [1, 1, 1, 1, 1]), select=0):

    """
    factor: scale factor for the energy. useful for obtaining the local anisotropy of the 1st or 3rd Co atom.
    p0=(D, E, E0, dtheta, dphi)
    """

    ### Read data
    
    fin = base_name + ".dat"

    if not os.path.isfile(fin):
        return
    
    energies = np.loadtxt(fin)
    nconfig = energies.shape[0]



    ### Data massage

    columns = np.array(columns, dtype=int) - 1

    coli = columns[0]
    colj = columns[1]
    colk = columns[2]

    # Select data points for fitting.

    if select == 1:
        # Use the south and north poles and those on the equator.
        rows_to_remove = []
        for i in range(nconfig):
            if not ( energies[i, coli] == 0. or energies[i, coli] == 90. or energies[i, coli] == 180.):
                rows_to_remove.append(i)
        energies = np.delete(energies, rows_to_remove, 0)
        nconfig = energies.shape[0]
    else:
        # Use all data points
        pass

    emin = np.nanmin(energies[:, colk])

    energies = energies[np.argsort(energies[:,colk])]

    for i in range(nconfig):
        if energies[i, colk] == 0.0:
            energies[i, colk] = np.nan
        else:
            energies[i, colk] = 1000*factor*(energies[i, colk] - emin)
    
    emin = np.nanmin(energies[:, colk])
    emax = np.nanmax(energies[:, colk])

    energies[:, 0:colk] = np.deg2rad(energies[:, 0:colk])
    p0[3:5] = np.deg2rad(p0[3:5])
    bounds[0][3:5] = np.deg2rad(bounds[0][3:5])
    bounds[1][3:5] = np.deg2rad(bounds[1][3:5])
    
    

    ### Fit 

    popt, pcov = curve_fit(e_ani, (energies[:, coli], energies[:, colj]), energies[:, colk], p0=p0, bounds=bounds)
    print("D = {:12.6f} (meV), E = {:12.6f} (meV), E0 = {:12.6f} (meV), δθ = {:6.1f} (deg), δφ = {:6.1f} (deg)\n".format(*popt[0:3], *np.rad2deg(popt[3:5])))
    #print("D={:.3f} meV, E={:.3f} meV\n".format(*popt[0:2]))

    energies_fitted = e_ani((energies[:, coli], energies[:, colj]), D=popt[0], E=popt[1], E0=popt[2], dtheta=popt[3], dphi=popt[4])
    de = energies_fitted - energies[:, colk]
    energies[:, 0:colk] = np.rad2deg(energies[:, 0:colk])
    print("Maximal deviation in energy: {:12.6f} (meV)\n".format( np.max(np.abs(de)) ))
    print("Root mean square of energy differences: {:12.6f} (meV)\n".format( np.sqrt(np.mean(np.square(de))) ))

    with open("energies_fitted.dat", "w") as f:
        for i in range(nconfig):
            f.write("{:6.1f}  {:6.1f}  {:6.1f}   {:8.3f}  {:8.3f}\n".format(energies[i, coli], energies[i, colj], energies[i, colk], energies_fitted[i], de[i]))



    ### Local reference frame

    dxx = -popt[0]/3 + popt[1]
    dyy = -popt[0]/3 - popt[1]
    dzz = 2*popt[0]/3
    emat_loc = get_local_eigen_frame(popt[3], popt[4])

    # d_matrix v = w v, V = [[v11, v12, v13], [v21, v22, v23], [v31, v32, v33]]
    # V == emat_loc

    with open("d_spin-orbit.dat", "w") as f:
        f.write("{:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(dxx, *emat_loc[0]))
        f.write("{:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(dyy, *emat_loc[1]))
        f.write("{:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(dzz, *emat_loc[2]))



    ### Plot
    
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    
    mpl.rcParams["font.family"] = "Times New Roman"
    mpl.rcParams["font.size"] = 12
    
    ax.set_title(title, {'fontsize': 10})

    x = [i for i in range(nconfig)]
    
    ax.plot(x, energies[:,colk], 'bo', label="DFT")
    ax.plot(x, energies_fitted, 'rs', label="Fitted values")

    ax.legend()

    ax.set_xlabel("magnetic configuration")
    ax.set_ylabel("E (meV)")
    
    fig.savefig("fitting_quality" + ".pdf")
    

if __name__ == "__main__":

    # Reference frame in which angles are written for DFT energies.
    emat_dft = np.eye(3)

    e_ani = e_ani_q

    D_min      = -2.0     ;   D_max      = -0.2
    E_min      = -0.5     ;   E_max      = -0.001
    E0_min     = -10.0    ;   E0_max     =  10.0
    dtheta_min =  0.      ;   dtheta_max =  5.
    dphi_min   =  0.      ;   dphi_max   =  180.

    p0 = [ (D_min+D_max)/2, (E_min+E_max)/2, E0_min, dtheta_min, dphi_min+90]

    bounds=([D_min, E_min, E0_min, dtheta_min, dphi_min], [D_max, E_max, E0_max, dtheta_max, dphi_max])

    fit(base_name = "energies_average", columns=[1, 2, 3], title="Cr3", p0=p0, bounds=bounds, select=1)

