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
        data = image[1].data[image[3].data==0]
        data = data[np.logical_and(data<np.percentile(data,90),data>np.percentile(data,10))]
        skyavg[i] = np.average(data)
        skysig[i] = np.std(data)
        image.close()
    
    return skyavg, skysig