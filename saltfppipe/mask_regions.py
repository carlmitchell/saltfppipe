import numpy as np
from sys import exit
from saltfppipe.fp_image_class import FPImage


def mask_circle(array, par):

    # Make X and Y arrays
    xgrid, ygrid = np.meshgrid(np.arange(array.shape[1])+1,
                               np.arange(array.shape[0])+1)
    rad2grid = (xgrid-par[0])**2+(ygrid-par[1])**2
    array[rad2grid < par[2]**2] = 1

    return array


def mask_box(array, par):

    # Make X and Y arrays
    xgrid, ygrid = np.meshgrid(np.arange(array.shape[1])+1,
                               np.arange(array.shape[0])+1)

    # Get four vertices
    par[4] = np.mod(par[4], 180)
    if par[4] >= 90:
        par[4] -= 90
        par[2], par[3] = par[3], par[2]
    par[4] *= np.pi/180
    par[3] *= .5
    par[2] *= .5
    (v1x, v1y) = (par[0]+np.cos(par[4])*par[2]+np.sin(par[4])*par[3],
                  par[1]+np.sin(par[4])*par[2]-np.cos(par[4])*par[3])
    (v2x, v2y) = (par[0]+np.cos(par[4])*par[2]-np.sin(par[4])*par[3],
                  par[1]+np.sin(par[4])*par[2]+np.cos(par[4])*par[3])
    (v3x, v3y) = (par[0]-np.cos(par[4])*par[2]-np.sin(par[4])*par[3],
                  par[1]-np.sin(par[4])*par[2]+np.cos(par[4])*par[3])
    (v4x, v4y) = (par[0]-np.cos(par[4])*par[2]+np.sin(par[4])*par[3],
                  par[1]-np.sin(par[4])*par[2]-np.cos(par[4])*par[3])

    # Right boundary
    if v1x == v2x:
        mask = xgrid < v1x
    else:
        m = (v1y-v2y)/(v1x-v2x)
        b = v1y-m*v1x
        mask = xgrid < (ygrid-b)/m

    # Left boundary
    if v3x == v4x:
        mask = np.logical_and(mask, xgrid > v3x)
    else:
        m = (v3y-v4y)/(v3x-v4x)
        b = v3y-m*v3x
        mask = np.logical_and(mask, xgrid > (ygrid-b)/m)

    # Top boundary
    m = (v3y-v2y)/(v3x-v2x)
    b = v3y-m*v3x
    mask = np.logical_and(mask, ygrid < m*xgrid+b)

    # Bottom boundary
    m = (v1y-v4y)/(v1x-v4x)
    b = v1y-m*v1x
    mask = np.logical_and(mask, ygrid > m*xgrid+b)

    array[mask] = 1

    return array


def mask_ellipse(array, par):

    # Make coordinate arrays
    xgrid, ygrid = np.meshgrid(np.arange(array.shape[1])+1,
                               np.arange(array.shape[0])+1)
    anglegrid = np.arctan2(ygrid-par[1], xgrid-par[0])-par[4]*np.pi/180
    radgrid = np.sqrt((xgrid-par[0])**2+(ygrid-par[1])**2)
    newxgrid = radgrid*np.cos(anglegrid)
    newygrid = radgrid*np.sin(anglegrid)
    array[(newygrid/par[3])**2+(newxgrid/par[2])**2 < 1] = 1

    return array


def mask_regions(imagepath, regionpath):
    """Adds the ds9 regions from a specified region file to the bad pixel mask
    of an image.

    The region file should be of the kind output by ds9 (region->save regions)
    with the following modes selected:
    format: Anything other than 'XML' or 'X Y'
    coordinate system: 'image'

    Currently this code only supports the following region types:
    box, ellipse, circle

    Inputs:
    imagepath -> String, path to an image
    regionpath -> String, path to a region file

    """

    # Open and parse the region file
    infile = open(regionpath)
    a = infile.readlines()
    if "format: DS9" in a[0]:
        var = "ds9"
    elif "logical" in a[0]:
        var = "pros"
    elif "format: pixels" in a[0]:
        var = "tng"
    elif "format: CIAO" in a[0]:
        var = "ciao"
    elif "(" in a[0]:
        var = "saoimage"
    else:
        exit("Invalid region file.")
    shapes = []
    par = []
    if var == "ds9":
        for i in range(3, len(a)):
            line = a[i]
            line = line.replace(",", " ")
            line = line.replace("(", " ")
            line = line.replace(")", " ")
            line = line.split()
            shapes.append(line[0])
            for j in range(1, len(line)):
                line[j] = float(line[j])
            par.append(line[1:])
    if var == "pros":
        for i in range(len(a)):
            line = a[i].split()
            shapes.append(line[1])
            for j in range(2, len(line)):
                line[j] = float(line[j])
            par.append(line[2:])
    if var == "tng":
        for i in range(1, len(a)):
            line = a[i]
            line = line.replace("(", " ")
            line = line.replace(")", " ")
            line = line.replace("+", " ")
            line = line.replace(",", " ")
            line = line.split()
            shapes.append(line[0])
            for j in range(1, len(line)-2):
                line[j] = float(line[j])
            par.append(line[1:-2])
    if var == "ciao":
        for i in range(1, len(a)):
            line = a[i]
            line = line.replace("(", " ")
            line = line.replace(")", " ")
            line = line.replace(",", " ")
            line = line.replace("rotbox", "box")
            line = line.split()
            shapes.append(line[0])
            for j in range(1, len(line)):
                line[j] = float(line[j])
            par.append(line[1:])
    if var == "saoimage":
        for i in range(len(a)):
            line = a[i]
            line = line.replace("(", " ")
            line = line.replace(")", " ")
            line = line.replace(",", " ")
            line.split()
            shapes.append(line[0])
            for j in range(1, len(line)):
                line[j] = float(line[j])
            par.append(line[1:])

    # Open the input image and run it through the various region masks
    im = FPImage(imagepath, update=True)
    for i in range(len(shapes)):
        if shapes[i] == "circle":
            im.badp = mask_circle(im.badp, par[i])
        elif shapes[i] == "box":
            im.badp = mask_box(im.badp, par[i])
        elif shapes[i] == "ellipse":
            im.badp = mask_ellipse(im.badp, par[i])

    # Close files
    infile.close()
    im.close()

    return