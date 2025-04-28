import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def do(base_name = "energies", columns=(0, 1, 6), factor=1.0):
    # factor: scale factor for the energy. useful when multiple equivalent spins are rotated togethe.

    columns = np.array(columns, dtype=int) - 1

    ### Input and output
    
    fin = base_name + ".dat"
    fout = base_name + ".pdf"

    if not os.path.isfile(fin):
        return
    
    energies = np.loadtxt(fin)
    #print(energies[-5:]); exit()
    
    
    
    ### Data massage

    column1 = columns[0]
    column2 = columns[1]
    column3 = columns[2]
    
    emin = np.nanmin(energies[:, column3])
    #print("{:15.8f}".format(emin))
    
    nconfig = energies.shape[0]

    Ez = 0.0
    energies_theta90 = []
    for i in range(nconfig):
        if energies[i, column1] == 0.0 and energies[i, column2] == 0.0:
            Ez = energies[i, column3]

        if energies[i, column1] == 90.0:
            energies_theta90.append(energies[i, column3])
    
    easy_axis = False
    if Ez < energies_theta90[0]:
        easy_axis = True
    
    
    ### MAE

    if easy_axis:
        mae_out = 1000*(np.nanmax(energies_theta90) - Ez)
    else:
        mae_out = 1000*(Ez - np.nanmin(energies_theta90))
    mae_in  = 1000*(np.nanmax(energies_theta90) - np.nanmin(energies_theta90))
    
    print("Out-of-plane MAE: {:6.3f} meV.".format(mae_out))
    print("    In-plane MAE: {:6.3f} meV.".format(mae_in ))
    
if __name__ == "__main__":

    do(base_name = "energies_average", columns=(1,2,3), factor=1)

