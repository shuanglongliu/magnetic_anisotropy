import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl
import matplotlib.pyplot as plt

emat = np.eye(3)

def f(x, a, b, c):
    return a*(x-c)**2 + b

def fit(base_name = "energies.dat", columns=[1, 3], title="Parabola", p0=[0, 0, 0], bounds=([0, 0, 0], [1, 1, 1])):

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

    emin = np.nanmin(data[:, colj])

    for i in range(nconfig):
        if data[i, colj] == 0.0:
            data[i, colj] = np.nan
        else:
            data[i, colj] = 1000*(data[i, colj] - emin)
    
    
    ### Fit 

    popt, pcov = curve_fit(f, data[:, coli], data[:, colj], p0=p0, bounds=bounds)
    print("c = {:6.3f} (deg)\n".format(popt[2]))

    data_fitted = f(data[:, coli], a=popt[0], b=popt[1], c=popt[2])
    de = data_fitted - data[:, colj]
    print("Maximal deviation in energy: {:12.6f} (meV)\n".format( np.max(np.abs(de)) ))
    print("RMS deviation in energy: {:12.6f} (meV)\n".format( np.sqrt(np.mean(np.square(de))) ))


    ### Plot
    
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    
    mpl.rcParams["font.family"] = "Times New Roman"
    mpl.rcParams["font.size"] = 12
    
    ax.set_title(title, {'fontsize': 10})

    xmin = np.nanmin(data[:, coli])
    xmax = np.nanmax(data[:, coli])
    x = np.linspace(xmin, xmax, 100, endpoint=True)
    y = f(x, a=popt[0], b=popt[1], c=popt[2])

    ax.plot(data[:, coli], data[:,colj], 'bo', label="DFT")
    ax.plot(x, y, 'r-', label="Fitted values")

    ax.legend()

    ax.set_xlabel("θ (deg)")
    ax.set_ylabel("E (meV)")
    
    fig.tight_layout()
    fig.savefig("fitting_quality" + ".pdf")
    
if __name__ == "__main__":

    p0 = [1., 0., 1]

    bounds=([0., -5., -5.], [10., 0., 5.])

    fit(base_name = "energies", columns=[1, 3], title="", p0=p0, bounds=bounds)

