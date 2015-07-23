import numpy as np
from astropy.io.fits import open as openfits
from sys import exit as crash

def deghost(fn,uncfn=None,g=0.04):
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
    uncfn (optional) -> String, the path to the uncertainty image.
    g -> The multiple to subtract from the original image, default = 4%
    
    """
    
    #Open the image and check for center coordinates
    image = openfits(fn,mode="update")
    xcen = image[0].header.get("fpxcen")
    ycen = image[0].header.get("fpycen")
    if xcen == None:
        print "Error! Image "+fn+" doesn't have center coordinates in header!"
        crash()
    
    #Deghost the image
    print "Deghosting image "+fn
    image[0].header["fpghost"] = "True"
    
    #Determine image size
    xsize=image[0].data.shape[1]
    ysize=image[0].data.shape[0]

    #Make a mask for the chip gaps and outside-aperture stuff
    mask = image[0].data==0

    #Make an array of the flipped data
    flip = image[0].data[::-1,::-1].copy()

    #Calculate the difference between the image's geometric center (midpoint) and the axis of rotation
    xshift = 2*xcen-xsize-1
    yshift = 2*ycen-ysize-1

    #A given pixel's position, when rotated, will overlap its four neighboring pixels.
    #All pixels will have the same four overlapping regions because the rotation is 180 degrees.
    #Here we take a weighted sum of these four pixels, where the weights are equal to the areas of overlap.
    #This weighted sum is subtracted from the original pixel in the data array.
    image[0].data[max(np.floor(yshift),0):min(ysize+np.floor(yshift),ysize),max(np.floor(xshift),0):min(xsize+np.floor(xshift),xsize)] += ( -g*abs((np.ceil(yshift)-yshift)*(np.ceil(xshift)-xshift))*flip[max(-np.floor(yshift),0):min(ysize,ysize-np.floor(yshift)),max(-np.floor(xshift),0):min(xsize,xsize-np.floor(xshift))])
    image[0].data[max(np.ceil(yshift),0):min(ysize+np.ceil(yshift),ysize),max(np.floor(xshift),0):min(xsize+np.floor(xshift),xsize)] += ( -g*abs((np.floor(yshift)-yshift)*(np.ceil(xshift)-xshift))*flip[max(-np.ceil(yshift),0):min(ysize,ysize-np.ceil(yshift)),max(-np.floor(xshift),0):min(xsize,xsize-np.floor(xshift))])
    image[0].data[max(np.floor(yshift),0):min(ysize+np.floor(yshift),ysize),max(np.ceil(xshift),0):min(xsize+np.ceil(xshift),xsize)] += ( -g*abs((np.ceil(yshift)-yshift)*(np.floor(xshift)-xshift))*flip[max(-np.floor(yshift),0):min(ysize,ysize-np.floor(yshift)),max(-np.ceil(xshift),0):min(xsize,xsize-np.ceil(xshift))])
    image[0].data[max(np.ceil(yshift),0):min(ysize+np.ceil(yshift),ysize),max(np.ceil(xshift),0):min(xsize+np.ceil(xshift),xsize)] += ( -g*abs((np.floor(yshift)-yshift)*(np.floor(xshift)-xshift))*flip[max(-np.ceil(yshift),0):min(ysize,ysize-np.ceil(yshift)),max(-np.ceil(xshift),0):min(xsize,xsize-np.ceil(xshift))])
    
    #Remask the data using the mask
    image[0].data[mask] = 0
    
    #Close the image
    image.close()
    
    #Modify the uncertainty image if it exists
    if not (uncfn is None):
        print "Updating uncertainty image "+uncfn
        uncimage = openfits(uncfn,mode="update")
        
        #Make an array of the flipped data
        flip = uncimage[0].data[::-1,::-1].copy()*1.0
        
        #Uncertainties add in quadrature
        uncimage[0].data = np.power(uncimage[0].data,2)
        
        #Add the uncertainty in the deghosted uncertainty image to the original
        uncimage[0].data[max(np.floor(yshift),0):min(ysize+np.floor(yshift),ysize),max(np.floor(xshift),0):min(xsize+np.floor(xshift),xsize)] += (g*abs((np.ceil(yshift)-yshift)*(np.ceil(xshift)-xshift))*flip[max(-np.floor(yshift),0):min(ysize,ysize-np.floor(yshift)),max(-np.floor(xshift),0):min(xsize,xsize-np.floor(xshift))])**2
        uncimage[0].data[max(np.ceil(yshift),0):min(ysize+np.ceil(yshift),ysize),max(np.floor(xshift),0):min(xsize+np.floor(xshift),xsize)] += ( g*abs((np.floor(yshift)-yshift)*(np.ceil(xshift)-xshift))*flip[max(-np.ceil(yshift),0):min(ysize,ysize-np.ceil(yshift)),max(-np.floor(xshift),0):min(xsize,xsize-np.floor(xshift))])**2
        uncimage[0].data[max(np.floor(yshift),0):min(ysize+np.floor(yshift),ysize),max(np.ceil(xshift),0):min(xsize+np.ceil(xshift),xsize)] += ( g*abs((np.ceil(yshift)-yshift)*(np.floor(xshift)-xshift))*flip[max(-np.floor(yshift),0):min(ysize,ysize-np.floor(yshift)),max(-np.ceil(xshift),0):min(xsize,xsize-np.ceil(xshift))])**2
        uncimage[0].data[max(np.ceil(yshift),0):min(ysize+np.ceil(yshift),ysize),max(np.ceil(xshift),0):min(xsize+np.ceil(xshift),xsize)] += ( g*abs((np.floor(yshift)-yshift)*(np.floor(xshift)-xshift))*flip[max(-np.ceil(yshift),0):min(ysize,ysize-np.ceil(yshift)),max(-np.ceil(xshift),0):min(xsize,xsize-np.ceil(xshift))])**2
        
        uncimage[0].data = np.sqrt(uncimage[0].data)
        
        #Close the uncertainty image
        uncimage.close()
        
    return