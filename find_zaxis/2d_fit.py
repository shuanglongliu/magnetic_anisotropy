import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.spatial.transform import Rotation as R
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylib import change_frame_sph_rad

emat_dft = np.eye(3)

def rotate_frame(emat=np.eye(3), dtheta=0., dphi=0.):

    """
    Rotate the local reference frame
      rotate the frame around the y axis by detheta, and then
      rotate the the frame around the old z axis by dphi
    Unit for angle: rad
    """

    # =================================================
    # rotate the frame around the y axis by detheta
    # =================================================

    r = R.from_rotvec( dtheta * emat[1] )
    rotmat = r.as_matrix()
    ex = np.matmul(rotmat, emat[0])
    ey = emat[1]
    ez = np.matmul(rotmat, emat[2])

    # =======================================================
    # rotate the frame around the old z axis by dphi
    # =======================================================

    r = R.from_rotvec( dphi * emat[2] )
    rotmat = r.as_matrix()
    ex = np.matmul(rotmat, ex)
    ey = np.matmul(rotmat, ey)
    ez = np.matmul(rotmat, ez)

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


def f(angles, a, c, dtheta, dphi):

    global emat_dft

    thetas, phis = change_frame(angles, emat_dft, dtheta, dphi)

    return a*thetas**2 + c

def fit(base_name = "energies.dat", columns=[1, 2, 3], title="Parabola", p0=[0, 0, 0], bounds=([0, 0, 0], [1, 1, 1])):

    ### Read data
    
    fin = base_name + ".dat"

    if not os.path.isfile(fin):
        return
    
    data = np.loadtxt(fin)
    nconfig = data.shape[0]



    ### Data massage
    
    columns = np.array(columns, dtype=int) - 1

    coli = columns[0]
    colj = columns[1]
    colk = columns[2]

    emin = np.nanmin(data[:, colk])

    data = data[np.argsort(data[:,colk])]

    for i in range(nconfig):
        if data[i, colk] == 0.0:
            data[i, colk] = np.nan
        else:
            data[i, colk] = 1000*(data[i, colk] - emin)
    
    data[:, 0:colk] = np.deg2rad(data[:, 0:colk])
    p0[2:4] = np.deg2rad(p0[2:4])
    bounds[0][2:4] = np.deg2rad(bounds[0][2:4])
    bounds[1][2:4] = np.deg2rad(bounds[1][2:4])
    
    
    ### Fit 

    popt, pcov = curve_fit(f, (data[:, coli], data[:, colj]), data[:, colk], p0=p0, bounds=bounds)
    print("dtheta = {:6.3f} (deg), dphi = {:6.3f} (deg)\n".format(np.rad2deg(popt[2]), np.rad2deg(popt[3])))

    data_fitted = f((data[:, coli], data[:, colj]), a=popt[0], c=popt[1], dtheta=popt[2], dphi=popt[3])
    de = data_fitted - data[:, colk]
    print("Maximal deviation in energy: {:12.6f} (meV)\n".format( np.max(np.abs(de)) ))
    print("RMS deviation in energy: {:12.6f} (meV)\n".format( np.sqrt(np.mean(np.square(de))) ))


    ### Plot
    
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    
    mpl.rcParams["font.family"] = "Times New Roman"
    mpl.rcParams["font.size"] = 12
    
    ax.set_title(title, {'fontsize': 10})

    xmin = np.nanmin(data[:, coli])
    xmax = np.nanmax(data[:, coli])

    x = [i for i in range(nconfig)]

    ax.plot(x, data[:,colk], 'bo', label="DFT")
    ax.plot(x, data_fitted[:], 'rs', label="Fitted values")

    ax.legend()

    ax.set_ylabel("E (meV)")
    
    fig.tight_layout()
    fig.savefig("fitting_quality" + ".pdf")
    
if __name__ == "__main__":

    p0 = [1., 0., 1, 45]

    bounds=([0., -5., -5., 0.], [10., 5., 5., 360.])

    fit(base_name = "energies", columns=[1, 2, 3], title="CrFe2 SMM SCO", p0=p0, bounds=bounds)

