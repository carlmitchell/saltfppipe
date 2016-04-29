from astropy.io import fits
import numpy as np
from saltfppipe.fp_image_class import FPImage


def combine_flat(flatlist, outfilepath="flat.fits"):
    """Combines a series of flatfield exposures to create a single flatfield
    image. Uses only the information contained inside the RSS aperture, which
    is the only significant difference between this routine and the familiar
    flat combination routines common to CCD astronomy.

    Inputs:
    flatlist -> A list of strings, each the path to a flatfield exposure
    outfilepath -> A string, the path of the output file (default 'flat.fits')

    """

    images = []
    for i in range(len(flatlist)):
        images.append(FPImage(flatlist[i]))
        rgrid = images[i].rarray(images[i].axcen,
                                 images[i].aycen)
        images[i].inty /= np.median(images[i].inty[rgrid < images[i].arad])
    imagestack = np.empty((len(flatlist), images[0].ysize, images[0].xsize))
    for i in range(len(flatlist)):
        imagestack[i, :, :] = images[i].inty
    flatarray = np.median(imagestack, axis=0)
    flatarray[flatarray < 0.000001] = 0
    fits.writeto(outfilepath, flatarray, clobber=True)

    # Close images
    for i in range(len(flatlist)):
        images[i].close()

    return
