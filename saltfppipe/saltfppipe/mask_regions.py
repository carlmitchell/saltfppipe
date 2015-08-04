import numpy as np
from astropy.io.fits import open as openfits
from sys import exit as crash

def mask_circle(array,params):
    xgrid,ygrid = np.meshgrid(np.arange(array.shape[1])+1,np.arange(array.shape[0])+1)
    rad2grid = (xgrid-params[0])**2+(ygrid-params[1])**2
    array[rad2grid<params[2]**2] = 1
    
    return array

def mask_box(array,params):
    
    #Make X and Y arrays
    xgrid,ygrid = np.meshgrid(np.arange(array.shape[1])+1,np.arange(array.shape[0])+1)
    
    #Get four vertices
    params[4]=np.mod(params[4],180)
    if params[4]>=90:
        params[4]-=90
        params[2],params[3]=params[3],params[2]    
    params[4]*=np.pi/180
    params[3]*=.5
    params[2]*=.5
    (v1x, v1y) = (params[0]+np.cos(params[4])*params[2]+np.sin(params[4])*params[3],
                  params[1]+np.sin(params[4])*params[2]-np.cos(params[4])*params[3])
    (v2x, v2y) = (params[0]+np.cos(params[4])*params[2]-np.sin(params[4])*params[3],
                  params[1]+np.sin(params[4])*params[2]+np.cos(params[4])*params[3])
    (v3x, v3y) = (params[0]-np.cos(params[4])*params[2]-np.sin(params[4])*params[3],
                  params[1]-np.sin(params[4])*params[2]+np.cos(params[4])*params[3])
    (v4x, v4y) = (params[0]-np.cos(params[4])*params[2]+np.sin(params[4])*params[3],
                  params[1]-np.sin(params[4])*params[2]-np.cos(params[4])*params[3])
    
    #Right boundary
    if v1x==v2x: mask = xgrid<v1x
    else:
        m = (v1y-v2y)/(v1x-v2x)
        b = v1y-m*v1x
        mask = xgrid<(ygrid-b)/m
    
    #Left boundary
    if v3x==v4x: mask = np.logical_and(mask,xgrid>v3x)
    else:
        m = (v3y-v4y)/(v3x-v4x)
        b = v3y-m*v3x
        mask = np.logical_and(mask,xgrid>(ygrid-b)/m)
    
    #Top boundary
    m = (v3y-v2y)/(v3x-v2x)
    b = v3y-m*v3x
    mask = np.logical_and(mask,ygrid<m*xgrid+b)
    
    #Bottom boundary
    m = (v1y-v4y)/(v1x-v4x)
    b = v1y-m*v1x
    mask = np.logical_and(mask,ygrid>m*xgrid+b)
    
    array[mask]=1
    
    return array

def mask_ellipse(array,params):
    xgrid,ygrid = np.meshgrid(np.arange(array.shape[1])+1,np.arange(array.shape[0])+1)
    anglegrid = np.arctan2(ygrid-params[1],xgrid-params[0])-params[4]*np.pi/180
    radgrid = np.sqrt((xgrid-params[0])**2+(ygrid-params[1])**2)
    newxgrid = radgrid*np.cos(anglegrid)
    newygrid = radgrid*np.sin(anglegrid)
    array[(newygrid/params[3])**2+(newxgrid/params[2])**2<1]=1

    return array

def mask_regions(imagepath,regionpath):
    """Adds the ds9 regions from a specified region file to the bad pixel mask
    of an image.
    
    The region file should be of the kind output by ds9 (region -> save regions)
    with the following modes selected:
    format: Anything other than 'XML' or 'X Y'
    coordinate system: 'image'
    
    Currently this code only supports the following region types:
    box, ellipse, circle
    
    Inputs:
    imagepath -> String, path to an image
    regionpath -> String, path to a region file
    
    """
    
    #Open and parse the region file
    infile = open(regionpath)
    a = infile.readlines()
    if "format: DS9" in a[0]: var = "ds9"
    elif "logical" in a[0]: var = "pros"
    elif "format: pixels" in a[0]: var = "tng"
    elif "format: CIAO" in a[0]: var = "ciao"
    elif "(" in a[0]: var = "saoimage"
    else: crash("Invalid region file.")
    shapes = []
    params = []
    if var == "ds9":
        for i in range(3,len(a)):
            line = a[i].replace("("," ").replace(")"," ").replace(","," ").split()
            shapes.append(line[0])
            for j in range(1,len(line)): line[j]=float(line[j])
            params.append(line[1:])
    if var == "pros":
        for i in range(len(a)):
            line = a[i].split()
            shapes.append(line[1])
            for j in range(2,len(line)): line[j]=float(line[j])
            params.append(line[2:])
    if var == "tng":
        for i in range(1,len(a)):
            line = a[i].replace("("," ").replace(")"," ").replace("+"," ").replace(","," ").split()
            shapes.append(line[0])
            for j in range(1,len(line)-2): line[j]=float(line[j])
            params.append(line[1:-2])
    if var == "ciao":
        for i in range(1,len(a)):
            line = a[i].replace("("," ").replace(")"," ").replace(","," ").replace("rotbox","box").split()
            shapes.append(line[0])
            for j in range(1,len(line)): line[j]=float(line[j])
            params.append(line[1:])
    if var == "saoimage":
        for i in range(len(a)):
            line = a[i].replace("("," ").replace(")"," ").replace(","," ").split()
            shapes.append(line[0])
            for j in range(1,len(line)): line[j]=float(line[j])
            params.append(line[1:])
    
    #Open the input image and run it through the various region masks
    im = openfits(imagepath,mode="update")
    for i in range(len(shapes)):
        if shapes[i]=="circle": im[3].data = mask_circle(im[3].data,params[i])
        elif shapes[i]=="box": im[3].data = mask_box(im[3].data,params[i])
        elif shapes[i]=="ellipse": im[3].data = mask_ellipse(im[3].data,params[i])
    
    #Close files
    infile.close()
    im.close()
    
    return