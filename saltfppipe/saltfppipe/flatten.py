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
        image[1].data *= 1/flatimage[0].data
        image[3].data *= 1/flatimage[0].data**2
        image[3].data[np.isnan(image[1].data)] = 1
        image[1].data[np.isnan(image[1].data)] = 0
        image[3].data[np.isinf(image[1].data)] = 1
        image[1].data[np.isinf(image[1].data)] = 0
        image.close()
    flatimage.close()