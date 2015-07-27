from astropy.io.fits import open as openfits
from astropy.io.fits import writeto as writefits
import numpy as np

def combine_flat(flatlist,outfilepath="flat.fits"):
    """Combines a series of flatfield exposures to create a single flatfield
    image. Uses only the information contained inside the RSS aperture, which
    is the only significant difference between this routine and the familiar
    flat combination routines common to CCD astronomy.
    
    Inputs:
    flatlist -> A list of strings, each the path to a flatfield exposure
    outfilepath -> A string, the path of the output file (default 'flat.fits')
    
    """
    
    imagelist = []
    for i in range(len(flatlist)):
        imagelist.append(openfits(flatlist[i]))
        xgrid, ygrid = np.meshgrid(np.arange(imagelist[i][1].data.shape[1]),np.arange(imagelist[i][1].data.shape[0]))
        imagelist[i][1].data=imagelist[i][1].data/np.median(imagelist[i][1].data[np.power((xgrid-imagelist[i][0].header["fpaxcen"]),2)+np.power((ygrid-imagelist[i][0].header["fpaycen"]),2)<np.power(imagelist[i][0].header["fparad"],2)])
    imagestack = np.empty((len(flatlist),imagelist[0][1].data.shape[0],imagelist[0][1].data.shape[1]))
    for i in range(len(flatlist)):
        imagestack[i,:,:] = imagelist[i][1].data
    flatarray = np.median(imagestack,axis=0)
    flatarray[flatarray<0.000001]=0
    writefits(outfilepath,flatarray,clobber=True)
    
    #Close images
    for i in range(len(flatlist)): imagelist[i].close()
    
    return