from saltfppipe.fp_image_class import FPImage


def aperture_mask(fnlist, axcen, aycen, arad):
    """Opens a series of images and masks pixels outside a circular aperture.
    The images are overwritten with the new masked image. The values of axcen,
    aycen, and arad are added to the FITS headers as keywords 'fpaxcen',
    'fpaycen', and 'fparad' respectively.

    Inputs:
    fnlist -> A list containing strings, each the path to a fits image.
    axcen -> The aperture x center
    aycen -> The aperture y center
    arad -> The aperture radius

    """

    for i in range(len(fnlist)):
        print ("Masking pixels outside aperture for image " +
               str(i+1)+" of "+str(len(fnlist))+": "+fnlist[i])
        image = FPImage(fnlist[i], update=True)
        rgrid = image.rarray(axcen, aycen)
        image.inty[rgrid > arad] = 0
        image.badp[rgrid > arad] = 1
        image.axcen = axcen
        image.aycen = aycen
        image.arad = arad
        image.close()
