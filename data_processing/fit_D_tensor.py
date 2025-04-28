import os
import sys
import numpy as np
from scipy.optimize import curve_fit
from scipy.spatial.transform import Rotation as R
import matplotlib as mpl
import matplotlib.pyplot as plt
from spin import Spin, dim, Sx, Sy, Sz, ID, ss, Sx2, Sy2, Sz2, SxSy, SySx, SxSz, SzSx, SySz, SzSy
from pylib import change_frame_sph_rad, sph2cart_rad, cart2sph_deg

def e_ani_q(angles, *params):

    """
    Magnetic anisotropy energy by treating spins as quantum operators.
    Spins are treated as operators.
    Units for angle: rad
    """

    thetas, phis = angles
    Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, E0 = params

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
    
        indices = np.argsort(w)
        i = dim-1
        state = v[:, indices[i]]

        op = Dxx*Sx2 + Dyy*Sy2 + Dzz*Sz2 + Dxy*(SxSy + SySx) + Dxz*(SxSz + SzSx) + Dyz*(SySz + SzSy)
        energy = np.dot(np.conjugate(state), np.matmul(op, state))

        energy = energy + E0

        energies.append(energy)
    energies = np.real(energies)

    return energies

def fit(base_name = "energies.dat", columns=[1, 2, 3], factor=1.0, title="mol", p0=[0, 0, 0, 0, 0, 0, 0], bounds=([0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1])):

    """
    factor: scale factor for the energy. 
    p0=(Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, E0)
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
    
    

    ### Fit 

    popt, pcov = curve_fit(e_ani, (energies[:, coli], energies[:, colj]), energies[:, colk], p0=p0, bounds=bounds)

    print()
    print("The fitted D tensor matrix elements are: ")
    print("Dxx = {:12.6f} (meV)".format(popt[0]))
    print("Dyy = {:12.6f} (meV)".format(popt[1]))
    print("Dzz = {:12.6f} (meV)".format(popt[2]))
    print("Dxy = {:12.6f} (meV)".format(popt[3]))
    print("Dxz = {:12.6f} (meV)".format(popt[4]))
    print("Dyz = {:12.6f} (meV)".format(popt[5]))
    print(" E0 = {:12.6f} (meV)".format(popt[6]))

    energies_fitted = e_ani((energies[:, coli], energies[:, colj]), *popt)
    de = energies_fitted - energies[:, colk]

    print()
    print("Maximal deviation in energy: {:12.6f} (meV).".format( np.max(np.abs(de)) ))
    print("Root mean square of energy differences: {:12.6f} (meV)".format( np.sqrt(np.mean(np.square(de))) ))

    energies[:, 0:colk] = np.rad2deg(energies[:, 0:colk])
    with open("energies_fitted.dat", "w") as f:
        for i in range(nconfig):
            f.write("{:6.1f}  {:6.1f}   {:8.3f}  {:8.3f}  {:8.3f}\n".format(energies[i, coli], energies[i, colj], energies[i, colk], energies_fitted[i], de[i]))



    ### Local reference frame

    Dxx, Dyy, Dzz, Dxy, Dxz, Dyz = popt[0:6]

    D = np.array([[Dxx, Dxy, Dxz], [Dxy, Dyy, Dyz], [Dxz, Dyz, Dzz]])
    v, w = np.linalg.eigh(D)

    indices = np.argsort(v)

    if (4*v[indices[1]] - 2*v[indices[2]] - 2*v[indices[0]] >= 0):
        # Easy plane
        i, j, k = (1, 2, 0)
    else:
        # Easy axis
        i, j, k = (1, 0, 2)

    Dxx = v[indices[i]]
    Dyy = v[indices[j]]
    Dzz = v[indices[k]]

    D = Dzz - (Dxx + Dyy)/2
    E = (Dxx - Dyy)/2

    print()
    print("In the local reference frame (eigenframe): ")
    print("D = {:9.3f} (meV), E = {:9.3f} (meV)".format(D, E))
    print("D = {:9.3f} (cm^-1), E = {:9.3f} (cm^-1)".format(8.06554*D, 8.06554*E))

    ex = w[:, i]
    ey = w[:, j]
    ez = w[:, k]

    if ( np.dot( np.cross(ex, ey), ez) < 0 ):
        ez = -ez

    print()
    print("Basis vectors of the local reference frame: ")
    print("ex = ( {:12.6f}, {:12.6f}, {:12.6f} )".format(*ex))
    print("ey = ( {:12.6f}, {:12.6f}, {:12.6f} )".format(*ey))
    print("ez = ( {:12.6f}, {:12.6f}, {:12.6f} )".format(*ez))

    print()



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

    Dxx_min      = -2.0     ;   Dxx_max      =  2.0
    Dyy_min      = -2.0     ;   Dyy_max      =  2.0
    Dzz_min      = -2.0     ;   Dzz_max      =  2.0
    Dxy_min      = -2.0     ;   Dxy_max      =  2.0
    Dxz_min      = -2.0     ;   Dxz_max      =  2.0
    Dyz_min      = -2.0     ;   Dyz_max      =  2.0
    E0_min     = -10.0    ;   E0_max     =  10.0

    p0 = [ (Dxx_min+Dxx_max)/2, (Dyy_min+Dyy_max)/2, (Dzz_min+Dzz_max)/2, (Dxy_min+Dxy_max)/2, (Dxz_min+Dxz_max)/2, (Dyz_min+Dyz_max)/2, E0_min]

    bounds=([Dxx_min, Dyy_min, Dzz_min, Dxy_min, Dxz_min, Dyz_min, E0_min], \
            [Dxx_max, Dyy_max, Dzz_max, Dxy_max, Dxz_max, Dyz_max, E0_max])

    fit(base_name = "energies_average", columns=[1, 2, 3], title="Pair1 of Mn8", p0=p0, bounds=bounds)

