from astropy.io.fits import open as openfits
import numpy as np

def flatten(fnlist,flatfile):
    """Flattens a series of images by dividing them by a flatfield image.
    This differs from traditional flatfielding only in that in maintains the
    RSS aperture mask.
    
    Inputs:
    fnlist -> A list of strings, each the path to an image to flatten
    flatfile -> The location of a master flat image
    
    """
    
    flatimage = openfits(flatfile)
    for i in range(len(fnlist)):
        print "Flattening image "+str(i+1)+" of "+str(len(fnlist))+": "+fnlist[i]
        image = openfits(fnlist[i],mode="update")
        image[0].header["fpflat"] = "True"
        xgrid, ygrid = np.meshgrid(np.arange(image[0].data.shape[1]),np.arange(image[0].data.shape[0]))
        image[0].data[np.power((xgrid-image[0].header["fpaxcen"]),2)+np.power((ygrid-image[0].header["fpaycen"]),2)<np.power(image[0].header["fparad"],2)] *= 1/flatimage[0].data[np.power((xgrid-image[0].header["fpaxcen"]),2)+np.power((ygrid-image[0].header["fpaycen"]),2)<np.power(image[0].header["fparad"],2)]
        image[0].data[np.isnan(image[0].data)] = 0
        image[0].data[np.isinf(image[0].data)] = 0
        image.close()
    flatimage.close()