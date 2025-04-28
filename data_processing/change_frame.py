import sys
import numpy as np

def sph2cart_deg(sph):
    # Convert cartesian coordinate to spherical coordinate
    sph[1:3] = np.deg2rad(sph[1:3])
    z = sph[0] * math.cos( sph[1] )
    x = sph[0] * math.sin( sph[1] ) * math.cos( sph[2] )
    y = sph[0] * math.sin( sph[1] ) * math.sin( sph[2] )
    return [x, y, z]

def cart2sph_deg(xyz):
    # Convert cartesian coordinate to spherical coordinate
    sph = []
    ptsnew = appendSpherical_np(np.array([xyz]))
    sph.append(ptsnew[0][3])
    sph.append(ptsnew[0][4]/math.pi*180)
    sph.append(ptsnew[0][5]/math.pi*180)
    return sph

def change_frame_sph_deg(vi_sph, emati, ematf):
    """
    emati: initial basis vectors 
    ematf: final basis vectors 
    vi_sph: the vector in the initial frame 
    vi_sph: the same vector in the final frame 
    """
    
    emati_t = np.transpose(emati)
    ematf_t = np.transpose(ematf)
    ematf_t_inv = np.linalg.inv(ematf_t)
    
    vi_cart = sph2cart_deg(vi_sph)
    vf_cart = np.matmul(emati_t, vi_cart)
    vf_cart = np.matmul(ematf_t_inv, vf_cart)
    
    vf_sph = cart2sph_deg(vf_cart)

    return vf_sph

def get_emat_local(z_direction=(0,0)):
    # let the local z axis be along the z_direction

    theta, phi = z_direction

    if abs(theta-90) < 1e-6 and abs(phi) < 1e-6:
        # local z axis is along the global x axis
        # swap axes: x -> z, y -> x, z -> y
        emat_local = [[0,0,1], [1,0,0], [0,1,0]]
    else:
        # define local x axis such that the global x axis is in the local xz plane
        ez = sph2cart_deg([1, z_direction[0], z_direction[1]])
        ey = np.cross(ez, [1, 0, 0])
        ey = ey/np.linalg.norm(ey)
        ex = np.cross(ey, ez)
        emat_local = [ex, ey, ez]
    emat_local = np.array(emat_local)

    return emat_local

def change_frame(base_name = "energies", columns=(2,3,6), z_direction=(0., 0.), z_direction_from_data=0):
    
    fin = base_name + ".dat"
    fout = "new_frame.dat"

    es = np.loadtxt(fin)
    ne = es.shape[0]
    nsite = int( (es.shape[1] - 1)/2 )
    #print(nsite, type(nsite)); exit()


    column1 = columns[0] - 1
    column2 = columns[1] - 1
    column3 = columns[2] - 1
    
    es = es[np.argsort(es[:,column3])]
    #print(es[0])

    if z_direction_from_data == 1:
        # Easy axis
        theta = es[0, column1]
        phi = es[0, column2]
        z_direction = (theta, phi)
    elif z_direction_from_data == 2:
        # Hard axis
        theta = es[-1, column1]
        phi = es[-1, column2]
        z_direction = (theta, phi)

    emat_new = get_emat_local(z_direction)

    es_new = []
    for i in range(ne):
        sph_old = [1, es[i,column1], es[i, column2]]
        r, theta, phi = change_frame_sph_deg(sph_old, np.eye(3), emat_new)
        es_new.append([theta, phi, es[i, column3]])
    es_new = np.array(es_new)
    #print(es_new[0])

    with open(fout, "w") as f:
        ostring = "{:6.1f} {:6.1f} {:15.9f}\n"
        for i in range(ne):
            f.write(ostring.format(*es_new[i]))

if __name__ == "__main__":

    # z_direction_from_data = 0: Use z_direction which is specified manually.
    # z_direction_from_data = 1: Use z_direction which corresponds to E_min.
    # z_direction_from_data = 2: Use z_direction which corresponds to E_max.

    change_frame(base_name = "energies_average", columns=(1,2,3), z_direction=(90, 60), z_direction_from_data=1)

