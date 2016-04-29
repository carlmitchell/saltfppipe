import numpy as np
from astropy.io.fits import open as openfits
from os.path import join


def make_clean_map(cubedir, clobber=False):
    """Applies a series of user-specified masks to a velocity map file. The
    current masks are:

    - Minimum S/N cutoff
    - Maximum velocity cutoff
    - Minimum velocity cutoff
    - Maximum chi^2 cutoff
    - Maximum dsigma/sigma cutoff
    - Maximum dinty/inty cutoff

    The values of these cutoffs are stored in the image headers of the output
    image.
    """

    # Open images
    velimage = openfits(join(cubedir, "velocity.fits"))
    intyimage = openfits(join(cubedir, "intensity.fits"))
    dintyimage = openfits(join(cubedir, "dintensity.fits"))
    dcontimage = openfits(join(cubedir, "dcontinuum.fits"))
    sigmaimage = openfits(join(cubedir, "sigma.fits"))
    dsigmaimage = openfits(join(cubedir, "dsigma.fits"))
    chi2image = openfits(join(cubedir, "chi2.fits"))

    # Derived arrays
    normalized_dsig = dsigmaimage[0].data/sigmaimage[0].data
    normalized_dinty = dintyimage[0].data/intyimage[0].data
    snarray = intyimage[0].data / dcontimage[0].data

    # User input for cuts to make
    sncut = float(raw_input("Enter minimum signal-to-noise: "))
    minvel = float(raw_input("Enter minimum velocity: "))
    maxvel = float(raw_input("Enter maximum velocity: "))
    maxchi2 = float(raw_input("Enter maximum chi^2 (or 0 for no cut): "))
    maxsigmarat = float(raw_input("Enter maximum dsigma/sigma ratio " +
                                  "(or 0 for no cut): "))
    maxintyrat = float(raw_input("Enter maximum dinty/inty ratio " +
                                 "(or 0 for no cut): "))

    # Mask for the cuts
    mask = snarray > sncut
    mask = np.logical_and(mask, velimage[0].data > minvel)
    mask = np.logical_and(mask, velimage[0].data < maxvel)
    mask = np.logical_and(mask, np.logical_not(np.isinf(snarray)))
    if maxchi2 != 0:
        mask = np.logical_and(mask, chi2image[0].data < maxchi2)
    if maxsigmarat != 0:
        mask = np.logical_and(mask, normalized_dsig < maxsigmarat)
    if maxintyrat != 0:
        mask = np.logical_and(mask, normalized_dinty < maxintyrat)

    # Sum of mask pixels
    print repr(np.sum(mask))+" pixels kept."

    # Mask the velocity map
    velimage[0].data[np.logical_not(mask)] = np.nan
    velimage[0].header["CUTSN"] = sncut
    velimage[0].header["CUTVEL1"] = minvel
    velimage[0].header["CUTVEL2"] = maxvel
    velimage[0].header["CUTCHI2"] = maxchi2
    velimage[0].header["CUTDSIG"] = maxsigmarat
    velimage[0].header["CUTDINTY"] = maxintyrat
    velimage.writeto(join(cubedir, "cleanvelmap.fits"), clobber=clobber)

    # Close images
    velimage.close()
    intyimage.close()
    dintyimage.close()
    dcontimage.close()
    sigmaimage.close()
    dsigmaimage.close()
    chi2image.close()

    return
