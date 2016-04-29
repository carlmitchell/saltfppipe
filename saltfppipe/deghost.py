import numpy as np
from sys import exit
from saltfppipe.fp_image_class import FPImage


def deghost(fn, g=0.04):
    """Routine to deghost an image by rotating it about a central point and
    subtracting a constant multiple of this rotated image from the original.
    My experience has shown that about 0.04 is the right value, but your
    mileage may vary.

    The fits header must contain values in "fpxcen" and "fpycen" for the center
    about which to be rotated.

    Optionally also modifies an uncertainty image to account for the effects
    of deghosting on the error propagation.

    Creates the fits header keyword "fpghost" = "True".

    Inputs:
    fn -> String, the path to the fits image to be deghosted.
    g -> The multiple to subtract from the original image, default = 4%

    """

    # Open the image and check for center coordinates
    image = FPImage(fn, update=True)
    if image.xcen is None:
        exit("Error! Image "+fn+" doesn't have center coordinates in header!")
    xcen = image.xcen+1
    ycen = image.ycen+1

    # Deghost the intensity image
    print "Deghosting image "+fn

    # Determine image size
    xsize, ysize = image.xsize, image.ysize

    # Make a mask for the chip gaps and outside-aperture stuff
    mask = image.inty == 0

    # Make an array of the flipped data
    flipinty = image.inty[::-1, ::-1].copy()*1.0
    flipvari = image.vari[::-1, ::-1].copy()*1.0

    # Calculate the difference between the image's geometric center (midpoint)
    # and the axis of rotation
    xshift = 2*xcen-xsize-1
    yshift = 2*ycen-ysize-1

    # A given pixel's position, when rotated, will overlap its four
    # neighboring pixels. All pixels will have the same four overlapping
    # regions because the rotation is 180 degrees. Here we take a weighted
    # sum of these four pixels, where the weights are equal to the areas
    # of overlap. This weighted sum is subtracted from the original pixel.
    for i in range(4):
        # This branching is for the four overlapping pixel regions
        if i == 0:
            miny = max(np.floor(yshift), 0)
            maxy = min(ysize+np.floor(yshift), ysize)
            minx = max(np.floor(xshift), 0)
            maxx = min(xsize+np.floor(xshift), xsize)
            flipminy = max(-np.floor(yshift), 0)
            flipmaxy = min(ysize, ysize-np.floor(yshift))
            flipminx = max(-np.floor(xshift), 0)
            flipmaxx = min(xsize, xsize-np.floor(xshift))
            frac = abs((np.ceil(yshift)-yshift)*(np.ceil(xshift)-xshift))
        if i == 1:
            miny = max(np.ceil(yshift), 0)
            maxy = min(ysize+np.ceil(yshift), ysize)
            minx = max(np.floor(xshift), 0)
            maxx = min(xsize+np.floor(xshift), xsize)
            flipminy = max(-np.ceil(yshift), 0)
            flipmaxy = min(ysize, ysize-np.ceil(yshift))
            flipminx = max(-np.floor(xshift), 0)
            flipmaxx = min(xsize, xsize-np.floor(xshift))
            frac = abs((np.floor(yshift)-yshift)*(np.ceil(xshift)-xshift))
        if i == 2:
            miny = max(np.floor(yshift), 0)
            maxy = min(ysize+np.floor(yshift), ysize)
            minx = max(np.ceil(xshift), 0)
            maxx = min(xsize+np.ceil(xshift), xsize)
            flipminy = max(-np.floor(yshift), 0)
            flipmaxy = min(ysize, ysize-np.floor(yshift))
            flipminx = max(-np.ceil(xshift), 0)
            flipmaxx = min(xsize, xsize-np.ceil(xshift))
            frac = abs((np.ceil(yshift)-yshift)*(np.floor(xshift)-xshift))
        if i == 3:
            miny = max(np.ceil(yshift), 0)
            maxy = min(ysize+np.ceil(yshift), ysize)
            minx = max(np.ceil(xshift), 0)
            maxx = min(xsize+np.ceil(xshift), xsize)
            flipminy = max(-np.ceil(yshift), 0)
            flipmaxy = min(ysize, ysize-np.ceil(yshift))
            flipminx = max(-np.ceil(xshift), 0)
            flipmaxx = min(xsize, xsize-np.ceil(xshift))
            frac = abs((np.floor(yshift)-yshift)*(np.floor(xshift)-xshift))

        # Rotate and subtract the intensity array
        image.inty[miny:maxy,
                   minx:maxx] -= g*frac*flipinty[flipminy:flipmaxy,
                                                 flipminx:flipmaxx]
        image.vari[miny:maxy,
                   minx:maxx] += (g*frac)**2*flipvari[flipminy:flipmaxy,
                                                      flipminx:flipmaxx]

    # Remask the intensity data using the mask we created
    image.inty[mask] = 0

    # Update the header
    image.ghosttog = "True"

    # Close the image
    image.close()

    return
