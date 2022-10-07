import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def polar_plot(base_name = "energies", title="", columns=(0, 1, 6), factor=1.0):
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
    
    for i in range(nconfig):
        if energies[i, column3] == 0.0:
            energies[i, column3] = np.nan
        else:
            energies[i, column3] = 1000*factor*(energies[i, column3] - emin)
    
    emin = np.nanmin(energies[:, column3])
    emax = np.nanmax(energies[:, column3])
    #print("{:10.5f} {:10.5f}".format(emin, emax)); exit()

    energies[:, column1] = np.deg2rad(energies[:, column1])
    energies[:, column2] = np.deg2rad(energies[:, column2])
    
    
    
    ### Plot
    
    fig = plt.figure(figsize=(5.5, 4.5))
    ax = fig.add_subplot(projection="polar")
    
    mpl.rcParams["font.family"] = "Times New Roman"
    mpl.rcParams["font.size"] = 12
    
    ax.set_title(title, {'fontsize': 13})
    
    ax.set_xlim([0, 2*np.pi])
    ax.set_ylim([0, np.pi])

    ax.set_yticks(np.linspace(0, np.pi, 7, endpoint=True), labels=[u"0\u00b0", u"30\u00b0", u"60\u00b0", u"90\u00b0", u"120\u00b0", u"150\u00b0", u"180\u00b0"])

    cm = mpl.cm.get_cmap('rainbow')
    
    scat = ax.scatter(energies[:,column2], energies[:,column1], s=60, c=energies[:,column3], vmin=emin, vmax=emax, cmap=cm)
    
    cbar = plt.colorbar(scat, label = u"E(\u03b8, \u03d5) (meV)", pad=0.1)

    fig.tight_layout()
    
    plt.savefig(fout)

if __name__ == "__main__":

    polar_plot(base_name = "energies", title="Cr3", columns=(1,2,3), factor=1)

