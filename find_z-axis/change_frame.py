import sys
import numpy as np
from pylib import sph2cart_deg, cart2sph_deg, change_frame_sph_deg, get_emat_local

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

