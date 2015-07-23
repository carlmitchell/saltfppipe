import numpy as np
from astropy.io.fits import open as openfits

def fit_sky_level(fnlist):
    """Fits for the average sky background and sky background noise in a series
    of images.
    
    Inputs:
    fnlist -> List of strings, each containing the path to a fits image.
    
    Outputs:
    skyavg -> Array containing the average sky level in each image
    skysig -> Array containing the st. dev. of the sky level in each image
    
    """
    
    skyavg = np.empty(len(fnlist))
    skysig = np.empty(len(fnlist))
    for i in range(len(fnlist)):
        image = openfits(fnlist[i])
        xgrid, ygrid = np.meshgrid(np.arange(image[0].data.shape[1]),np.arange(image[0].data.shape[0]))
        data = image[0].data[np.power((xgrid-image[0].header["fpaxcen"]),2)+np.power((ygrid-image[0].header["fpaycen"]),2)<np.power(image[0].header["fparad"],2)]
        data = data[np.logical_and(data<np.percentile(data,90),data>np.percentile(data,10))]
        skyavg[i] = np.average(data)
        skysig[i] = np.std(data)
        image.close()
    
    return skyavg, skysig