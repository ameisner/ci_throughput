import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

def plot_total_throughput(ps1=False):
    fname = '../etc/gfa_throughput-airmass_1.00.fits'
    tab = fits.getdata(fname)

    plt.plot(tab['LAMBDA_NM'], tab['THROUGHPUT'], c='k', linewidth=2)

    if ps1:
        fname_ps1 = '../etc/PAN-STARRS_PS1.r.fits'
        tab_ps1 = fits.getdata(fname_ps1)
        plt.plot(tab_ps1['LAMBDA_NM'], tab_ps1['THROUGHPUT'], c='r')
        
    plt.xlabel('wavelength (nm)', fontsize=14)
    plt.ylabel('total GFA throughput', fontsize=14)

    plt.xlim((np.min(tab['LAMBDA_NM']), np.max(tab['LAMBDA_NM'])))

    if not ps1:
        ymax = np.max(tab['THROUGHPUT']) + 0.03
    else:
        ymax = np.max(list(tab['THROUGHPUT']) + list(tab_ps1['THROUGHPUT'])) + 0.03
        plt.title('black = DESI GFA, red = PS1 r')
        
    plt.ylim((0, ymax))
    
    ax = plt.gca()
    ax.tick_params(axis="both", labelsize=14)

    ax.set_aspect(150.0)
    outname = 'gfa_total_throughput.eps'
    if ps1:
        outname = outname.replace('.eps', '-with_ps1.eps')

    plt.savefig(outname, bbox_inches='tight')
