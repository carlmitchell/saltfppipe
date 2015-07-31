import numpy as np
from astropy.io.fits import open as openfits
from astropy.io.fits import writeto as writefits
from os.path import join

def make_clean_map(cubedir, clobber=False):
    """
    """
    
    velimage = openfits(join(cubedir,"velocity.fits"))
    intyimage = openfits(join(cubedir,"intensity.fits"))
    dcontimage = openfits(join(cubedir,"dcontinuum.fits"))
    snarray = intyimage[0].data / dcontimage[0].data
    sncut = float(raw_input("Enter signal-to-noise cutoff: "))
    minvel = float(raw_input("Enter minimum velocity: "))
    maxvel = float(raw_input("Enter maximum velocity: "))
    mask = np.logical_and(snarray>sncut,np.logical_and(velimage[0].data>minvel, velimage[0].data<maxvel))
    print np.sum(mask)+" pixels kept."
    newvelarray = np.zeros_like(velimage[0].data)
    newvelarray[mask] = velimage[0].data[mask]
    writefits(join(cubedir,"cleanvelmap.fits"), newvelarray, clobber=clobber)
    velimage.close()
    intyimage.close()
    dcontimage.close()