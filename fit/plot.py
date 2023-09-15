import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def polar_plot(base_name = "energies", title="", columns=(0, 1, 6), factor=1.0, wavenumber=True):
    # factor: scale factor for the energy. useful when multiple equivalent spins are rotated together.

    columns = np.array(columns, dtype=int) - 1

    ### Input and output
    
    fin = base_name + ".dat"
    fout = base_name + ".pdf"

    if not os.path.isfile(fin):
        return
    
    energies = np.loadtxt(fin)
    
    
    
    ### Data massage

    column1 = columns[0]
    column2 = columns[1]
    column3 = columns[2]
    
    emin = np.nanmin(energies[:, column3])
    
    nconfig = energies.shape[0]
    
    for i in range(nconfig):
        if energies[i, column3] == 0.0:
            energies[i, column3] = np.nan
        else:
            if wavenumber:
                energies[i, column3] = 8.06554393738*1000*factor*(energies[i, column3] - emin)
            else:
                energies[i, column3] =               1000*factor*(energies[i, column3] - emin)
    
    emin = np.nanmin(energies[:, column3])
    emax = np.nanmax(energies[:, column3])

    energies[:, column1] = np.deg2rad(energies[:, column1])
    energies[:, column2] = np.deg2rad(energies[:, column2])

    
    
    ### Plot
    
    fig = plt.figure(figsize=(4.2, 3.2))
    ax = fig.add_subplot(projection="polar")
    
    mpl.rcParams["font.family"] = "Times New Roman"
    mpl.rcParams["font.size"] = 12
    
    ax.set_title(title, {'fontsize': 13})
    
    ax.set_xlim([0, 2*np.pi])
    ax.set_ylim([0, np.pi])

    ax.set_xticks(np.linspace(0, 2*np.pi, 9, endpoint=True), labels=[r"$\phi=0$", r"", r"$\phi=\pi/2$", r"", r"$\phi=\pi$", r"", r"$\phi=3\pi/2$", r"", r""])
    ax.set_yticks(np.linspace(0,   np.pi, 7, endpoint=True), labels=[r"$\theta=0$", r"", r"", r"$\theta=\pi/2$", r"", r"", r"$\theta=\pi$"])

    cm = mpl.cm.get_cmap('rainbow')
    
    scat = ax.scatter(energies[:,column2], energies[:,column1], s=60, c=energies[:,column3], vmin=emin, vmax=emax, cmap=cm)
    
    cbar = plt.colorbar(scat, label = r"$E(\theta, \phi)\; (\mathrm{cm}^{-1})$", pad=0.15)

    fig.tight_layout()
    
    plt.savefig(fout)

if __name__ == "__main__":

    polar_plot(base_name = "energies", title="", columns=(1,2,3), factor=1, wavenumber=True)

