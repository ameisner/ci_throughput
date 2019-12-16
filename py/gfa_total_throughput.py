import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

def plot_total_throughput():
    fname = '../etc/gfa_throughput-airmass_1.00.fits'
    tab = fits.getdata(fname)

    plt.plot(tab['LAMBDA_NM'], tab['THROUGHPUT'], c='k', linewidth=2)

    plt.xlabel('wavelength (nm)', fontsize=14)
    plt.ylabel('total GFA throughput', fontsize=14)

    plt.xlim((np.min(tab['LAMBDA_NM']), np.max(tab['LAMBDA_NM'])))

    plt.ylim((0, np.max(tab['throughput']) + 0.03))

    ax = plt.gca()
    ax.tick_params(axis="both", labelsize=14)

    ax.set_aspect(150.0)
    
    plt.savefig('gfa_total_throughput.eps', bbox_inches='tight')
