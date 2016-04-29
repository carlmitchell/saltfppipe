from pyraf import iraf
import sys
from cStringIO import StringIO
import numpy as np
from saltfppipe.fp_image_class import FPImage


def solar_velocity_shift(imagelist, restwave):
    """Opens the images of a data cube and shifts the wavelength planes
    so they are shifted to the solar frame. This is accomplished using the
    IRAF task 'rvcorrect'.

    It does it in a really really stupid way, which is to temporarily redirect
    IRAF's output to a string, because IRAF doesn't actually save the output
    anywhere (that I know of), it only prints it to the terminal. Then this
    routine searches through that string for the necessary value. It's a bit
    backwards, but it works!

    Note that the boost here is a Lorentz boost, not a Galilean boost, though
    values for either will be similar (i.e. the relativistic formula is used)

    Inputs:
    imagelist -> List of strings, the paths to wavelength images to be
                     updated.
    restwave -> The rest wavelength of the line you're interested in.

    """

    c = 299792.458  # in km/s

    for i in range(len(imagelist)):
        sys.stdout = iraf_output = StringIO()
        image = FPImage(imagelist[i], update=True)
        iraf.rvcorrect(header="N",
                       input="N",
                       imupdate="N",
                       observatory="SAAO",
                       year=float(image.datestring.split("-")[0]),
                       month=float(image.datestring.split("-")[1]),
                       day=float(image.datestring.split("-")[2]),
                       ut=image.ut,
                       ra=image.ra,
                       dec=image.dec,
                       vobs=0,
                       mode="a")
        image.solarvel = float(iraf_output.getvalue().split()[-6])
        iraf_output.close()
        beta_earth = (((image.wave/restwave)**2-1) /
                      ((image.wave/restwave)**2+1))
        beta_shift = image.solarvel/c
        beta_helio = (beta_earth+beta_shift)/(1+beta_earth*beta_shift)
        image.wave = restwave*(1+beta_helio)/np.sqrt(1-beta_helio**2)
        image.wave[np.isnan(image.wave)] = 0
        image.close()

    sys.stdout = sys.__stdout__

    return
