from astropy.io import fits
import numpy as np
from saltfppipe.fp_image_class import FPImage

def flatten(fnlist,flatfile):
    """Flattens a series of images by dividing them by a flatfield image.
    This differs from traditional flatfielding only in that in maintains the
    RSS aperture mask.
    
    Inputs:
    fnlist -> A list of strings, each the path to an image to flatten
    flatfile -> The location of a master flat image
    
    """
    
    flatimage = fits.open(flatfile)
    for i in range(len(fnlist)):
        print ( "Flattening image "+str(i+1)+
                " of "+str(len(fnlist))+": "+fnlist[i] )
        image = FPImage(fnlist[i],update=True)
        image.flattog = "True"
        image.inty /= flatimage[0].data       #Flatten data
        image.vari /= flatimage[0].data**2    #Flatten variance
        image.badp[np.isnan(image.inty)] = 1  #Fix bad pixels:
        image.inty[np.isnan(image.inty)] = 0
        image.badp[np.isinf(image.inty)] = 1
        image.inty[np.isinf(image.inty)] = 0
        image.close()
    flatimage.close()