import os
import sys
import numpy as np
from pylib import change_frame_sph, cart2sph, get_opposite_direction_sph

def get_energies_without_exchange(base_name = "energies", columns=[1, 2, 7]):
    # Take the average of E(direction) and E(opposite direction).

    columns = np.array(columns, dtype=int) - 1
    coli = columns[0]; colj = columns[1]; colk = columns[2]



    ### Input and output

    fin = base_name + ".dat"

    if not os.path.isfile(fin):
        return

    energies = np.loadtxt(fin)
    ne = len(energies)



    ### Convert the list of energies into a dictionary

    energies_dict = dict({})

    for i in range(ne):
        theta = "{:6.1f}".format(energies[i][coli])
        phi = "{:6.1f}".format(energies[i][colj])
        energies_dict[(theta, phi)] = energies[i][colk]

    #print(energies_dict)



    ### Take the average

    fout = base_name + "_average.dat"
    with open(fout, "w") as f:
        for i in range(ne):
            key1 = ("{:6.1f}".format(energies[i][coli]), "{:6.1f}".format(energies[i][colj]))
            theta, phi = get_opposite_direction_sph((energies[i][coli], energies[i][colj]))
            key2 = ("{:6.1f}".format(theta), "{:6.1f}".format(phi))
        
            if key1 in energies_dict.keys() and key2 in energies_dict.keys():
                e = (energies_dict[key1] + energies_dict[key2])/2
                f.write("{:9s} {:9s} {:20.12f}\n".format(key1[0], key1[1], e))
                #f.write("{:9s} {:9s} {:20.12f}\n".format(key2[0], key2[1], e))

if __name__ == "__main__":

    get_energies_without_exchange(base_name = "energies", columns=[1, 2, 7])

