import numpy as np
from astropy.io.fits import open as openfits
from astropy.io.fits import writeto as writefits
from sys import exit as crash
from saltfppipe.fit_sky_level import fit_sky_level

def convolve_uncert(uncertarry, intyarray, kern):
    """Convolves an uncertainty array by the kernel 'kern'.
    Ignores pixels with intensity values equal to zero so as to preserve
    user-done masking, chip gaps, aperture masking, etc.
    
    Note that the kernel must be square in shape, and furthermore the square
    must be odd in side length, i.e. 1x1, 3x3, 5x5, etc.
    
    Inputs:
    uncertarray -> Input wavelength array.
    intyarray -> Input intensity array. Masked values ought to be zero.
    kern -> Kernel for the convolution
    
    Outputs:
    final_uncert -> The convolved uncertainty array
    
    """
    
    #A mask for pixels which are nonzero
    mask = intyarray != 0
    
    #Kernel size for looping in the most straightforward way
    ksize = (kern.shape[1]-1)/2
    
    #Image dimensions
    xsize = intyarray.shape[1]
    ysize = intyarray.shape[0]
    
    #Empty uncertainty and weight arrays to fill while convolving
    running_unc_sum = np.zeros_like(intyarray)
    running_wgt_sum = np.zeros_like(intyarray)
    
    #Loop over the kernel
    for y in range(2*ksize+1):
        for x in range(2*ksize+1):
            running_unc_sum[ksize:ysize-ksize-1,ksize:xsize-ksize-1] += kern[y,x]**2*uncertarry[y:ysize-2*ksize-1+y,x:xsize-2*ksize-1+x]**2
            running_wgt_sum[ksize:ysize-ksize-1,ksize:xsize-ksize-1] += kern[y,x]*mask[y:ysize-2*ksize-1+y,x:xsize-2*ksize-1+x]
    
    #Divide the running uncertainty by the running kernel sum to normalize by
    #how much of the kernel we ended up using
    running_unc_sum[mask] = np.sqrt(running_unc_sum[mask])/running_wgt_sum[mask]
    running_unc_sum[np.logical_not(mask)]=0
    
    return running_unc_sum
    
def convolve_wave(wavearray, intyarray, kern):
    """Convolves a wavelength array by the kernel 'kern', weighting the
    convolution by intensity. Ignores pixels with intensity values equal to
    zero so as to preserve user-done masking, chip gaps, aperture masking, etc.
    
    Note that the kernel must be square in shape, and furthermore the square
    must be odd in side length, i.e. 1x1, 3x3, 5x5, etc.
    
    Inputs:
    wavearray -> Input wavelength array.
    intyarray -> Input intensity array. Masked values ought to be zero.
    kern -> Kernel for the convolution
    
    Outputs:
    final_wave -> The convolved wavelength array
    
    """    

    #A mask for pixels which are nonzero
    mask = intyarray != 0
    
    #Kernel size for looping in the most straightforward way
    ksize = (kern.shape[1]-1)/2
    
    #Image dimensions
    xsize = intyarray.shape[1]
    ysize = intyarray.shape[0]
    
    #Empty wavelength and weight arrays to fill while convolving
    running_wave_sum = np.zeros_like(intyarray)
    running_wgt_sum = np.zeros_like(intyarray)
    
    #Loop over the kernel
    for y in range(2*ksize+1):
        for x in range(2*ksize+1):
            running_wave_sum[ksize:ysize-ksize-1,ksize:xsize-ksize-1] += kern[y,x]*intyarray[y:ysize-2*ksize-1+y,x:xsize-2*ksize-1+x]*wavearray[y:ysize-2*ksize-1+y,x:xsize-2*ksize-1+x]
            running_wgt_sum[ksize:ysize-ksize-1,ksize:xsize-ksize-1] += kern[y,x]*intyarray[y:ysize-2*ksize-1+y,x:xsize-2*ksize-1+x]
    
    #Divide the running wavelength by the running weight sum to normalize by
    #how much of the kernel we ended up using
    running_wave_sum[mask] = running_wave_sum[mask]/running_wgt_sum[mask]
    running_wave_sum[np.logical_not(mask)]=0
    
    return running_wave_sum

def convolve_inty(intyarray,kern):
    """Convolves an intensity array by the kernel 'kern'. Ignores pixels with
    values equal to zero so as to preserve user-done masking, chip gaps,
    aperture masking, etc.
    
    Note that the kernel must be square in shape, and furthermore the square
    must be odd in side length, i.e. 1x1, 3x3, 5x5, etc.
    
    Inputs:
    intyarray -> Input intensity array. Masked values ought to be zero.
    kern -> Kernel for the convolution
    
    Outputs:
    final_inty -> The convolved intensity array
    
    """
    
    #A mask for pixels which are nonzero
    mask = intyarray != 0
    
    #Kernel size for looping in the most straightforward way
    ksize = (kern.shape[1]-1)/2
    
    #Image dimensions
    xsize = intyarray.shape[1]
    ysize = intyarray.shape[0]
    
    #Empty intensity and count arrays to fill while convolving
    running_inty_sum = np.zeros_like(intyarray)
    running_kern_sum = np.zeros_like(intyarray)
    
    #Loop over the kernel
    for y in range(2*ksize+1):
        for x in range(2*ksize+1):
            running_inty_sum[ksize:ysize-ksize-1,ksize:xsize-ksize-1] += kern[y,x]*intyarray[y:ysize-2*ksize-1+y,x:xsize-2*ksize-1+x]
            running_kern_sum[ksize:ysize-ksize-1,ksize:xsize-ksize-1] += kern[y,x]*mask[y:ysize-2*ksize-1+y,x:xsize-2*ksize-1+x]
    
    #Divide the running intensity by the running kernel sum to normalize by
    #how much of the kernel we ended up using
    running_inty_sum[mask] = running_inty_sum[mask]/running_kern_sum[mask]
    running_inty_sum[np.logical_not(mask)]=0
    
    return running_inty_sum

def make_final_image(input_image, output_image, output_wave_image,
                     desired_fwhm,
                     input_uncert_image=None, output_uncert_image=None,
                     clobber=False):
    """This routine makes the 'final' images for a data cube. At least the paths
    to the input image, output image, and output wavelength image are necessary
    for this. Beyond that, the user may also have the routine create uncertainty
    images as well.
    
    Images are convolved to the resolution 'desired_fwhm'. If the current fwhm
    is already higher than that, the routine will throw an error.
    
    A number of fits header keywords are necessary for this program to function
    properly. Any of these missing will throw an error.
    
    The output images are intensity-weighted, i.e. the wavelength image will be
    created such that the wavelengths at each pixel are the 'most likely'
    wavelength for the intensity at that pixel, etc.
    
    Inputs:
    input_image -> Path to the input image.
    output_image -> Path to the output image.
    output_wave_image -> Path to the output wavelength image.
    desired_fwhm -> Desired FWHM for the resultant image to have.
    
    Optional Inputs:
    input_uncert_image -> Path to the input uncertainty image, if it exists.
    output_uncert_image -> Path to the output uncertainty image, if it exists.
    clobber -> Overwrite output images if they already exist. Default is False.
    
    """
    
    print "Making final data cube images for image "+input_image
    
    #Measure the sky background level in the input image
    skyavg, skysig = fit_sky_level([input_image])
    
    #Open the input image and get various header keywords, crash if necessary
    intyimage = openfits(input_image)
    intygrid = intyimage[0].data
    fwhm = intyimage[0].header.get("fpfwhm")
    wave0 = intyimage[0].header.get("fpwave0")
    calf = intyimage[0].header.get("fpcalf")
    xcen = intyimage[0].header.get("fpxcen")
    ycen = intyimage[0].header.get("fpycen")
    if fwhm == None: crash("Error! FWHM not measured for image "+input_image+".")
    if wave0 == None or calf == None: crash("Error! Wavelength solution does "+
                                            "not exist for image "+input_image+".")
    if xcen == None or ycen == None: crash("Error! Center values not measured "+
                                           "image "+input_image+".")
    if fwhm>desired_fwhm: crash("Error! Desired FWHM too low for image "+
                                input_image+".")
    
    #Subtract the sky background from the image
    intygrid[intygrid!=0] -= skyavg[0]
    
    #Calculate the necessary FWHM for convolution and make the gaussian kernel
    fwhm_conv = np.sqrt(desired_fwhm**2-fwhm**2)
    sig = fwhm_conv/2.3548+0.0001
    ksize = np.ceil(4*sig) #Generate the kernel to 4-sigma
    kxgrid, kygrid = np.meshgrid(np.linspace(-ksize,ksize,2*ksize+1),np.linspace(-ksize,ksize,2*ksize+1))
    kern = np.exp(-(kxgrid**2+kygrid**2)/(2*sig**2)) #Gaussian (unnormalized because sig can be 0)
    kern = kern/np.sum(kern) #Normalize the kernel
    
    #Open and convolve the uncertainty image if one exists. Save the output.
    if input_uncert_image != None:
        uncertimage = openfits(input_uncert_image)
        uncertgrid = uncertimage[0].data
        #Add the sky background uncertainty to the uncertainty grid
        uncertgrid[intygrid!=0] = np.sqrt(uncertgrid[intygrid!=0]**2+skysig[0]**2)
        #Convolve the uncertainties appropriately
        new_uncert_grid = convolve_uncert(uncertgrid, intygrid, kern)
        #Write to output file
        writefits(output_uncert_image,new_uncert_grid,header=uncertimage[0].header,clobber=clobber)
        uncertimage.close()
        
    #Create and convolve the wavelength image. Save the output.
    xgrid, ygrid = np.meshgrid(np.arange(intyimage[0].data.shape[1]),
                               np.arange(intyimage[0].data.shape[0]))
    r2grid = (xgrid-xcen)**2 + (ygrid-ycen)**2
    wavegrid = wave0 / np.sqrt(1+r2grid/calf**2)
    newwavegrid = convolve_wave(wavegrid, intygrid, kern)
    writefits(output_wave_image,newwavegrid,header=intyimage[0].header,clobber=clobber)
    
    #Convolve the intensity image. Save the output
    newintygrid = convolve_inty(intygrid, kern)
    intyimage[0].header["fpfwhm"] = desired_fwhm #Update header FWHM keyword
    writefits(output_image,newintygrid,header=intyimage[0].header,clobber=clobber)
    
    #Close images
    intyimage.close()
    
    return