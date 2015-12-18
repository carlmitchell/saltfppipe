import numpy as np
from astropy.io.fits import ImageHDU
from sys import exit as crash
from saltfppipe.fp_image_class import FPImage

def convolve_variance(vararray, badpixarray, kern):
    """Convolves a variance array by the kernel 'kern'.
    Ignores pixels with intensity values equal to zero so as to preserve
    user-done masking, chip gaps, aperture masking, etc.
    
    Note that the kernel must be square in shape, and furthermore the square
    must be odd in side length, i.e. 1x1, 3x3, 5x5, etc.
    
    Inputs:
    vararray -> Input variance array.
    badpixarray -> Input bad pixel array. Bad values ought to be 1.
    kern -> Kernel for the convolution
    
    Outputs:
    final_variance -> The convolved uncertainty array
    
    """
    
    #A mask for good pixels
    mask = badpixarray!=1
    
    #Kernel size for looping in the most straightforward way
    ksize = (kern.shape[1]-1)/2
    
    #Image dimensions
    xsize = vararray.shape[1]
    ysize = vararray.shape[0]
    
    #Empty uncertainty and weight arrays to fill while convolving
    running_var_sum = np.zeros_like(vararray)
    running_wgt_sum = np.zeros_like(vararray)
    
    #Loop over the kernel
    for y in range(2*ksize+1):
        for x in range(2*ksize+1):
            ylim = ysize-2*ksize-1+y
            xlim = xsize-2*ksize-1+x
            running_var_sum[ksize:ysize-ksize-1,
                            ksize:xsize-ksize-1] += ( kern[y,x]**2 * 
                                                      vararray[y:ylim,x:xlim] * 
                                                      mask[y:ylim,x:xlim] )
            running_wgt_sum[ksize:ysize-ksize-1,
                            ksize:xsize-ksize-1] += ( kern[y,x]**2 * 
                                                      mask[y:ylim,x:xlim] )
    
    #Divide the running variance by the running weight sum to normalize by
    #how much of the kernel we ended up using
    running_var_sum[mask] = running_var_sum[mask]/running_wgt_sum[mask]
    running_var_sum[np.logical_not(mask)] = vararray[np.logical_not(mask)]
    
    return running_var_sum
    
def convolve_wave(wavearray, intyarray, badpixarray, kern):
    """Convolves a wavelength array by the kernel 'kern', weighting the
    convolution by intensity. Ignores pixels with intensity values equal to
    zero so as to preserve user-done masking, chip gaps, aperture masking, etc.
    
    Note that the kernel must be square in shape, and furthermore the square
    must be odd in side length, i.e. 1x1, 3x3, 5x5, etc.
    
    Inputs:
    wavearray -> Input wavelength array.
    intyarray -> Input intensity array. Masked values ought to be zero.
    badpixarray -> Bad pixel mask. Bad pixels ought to be 1.
    kern -> Kernel for the convolution
    
    Outputs:
    final_wave -> The convolved wavelength array
    
    """    
    
    #A mask for pixels which are good
    mask = np.logical_and(badpixarray==0,intyarray!=0)
    
    #Kernel size for looping in the most straightforward way
    ksize = (kern.shape[1]-1)/2
    
    #Image dimensions
    xsize = wavearray.shape[1]
    ysize = wavearray.shape[0]
    
    #Empty wavelength and weight arrays to fill while convolving
    running_wav_sum = np.zeros_like(wavearray)
    running_wgt_sum = np.zeros_like(wavearray)
    
    #Loop over the kernel
    for y in range(2*ksize+1):
        for x in range(2*ksize+1):
            ylim = ysize-2*ksize-1+y
            xlim = xsize-2*ksize-1+x
            running_wav_sum[ksize:ysize-ksize-1,
                            ksize:xsize-ksize-1] += ( kern[y,x] * 
                                                      intyarray[y:ylim,x:xlim] * 
                                                      wavearray[y:ylim,x:xlim] * 
                                                      mask[y:ylim,x:xlim] )
            running_wgt_sum[ksize:ysize-ksize-1,
                            ksize:xsize-ksize-1] += ( kern[y,x] * 
                                                      intyarray[y:ylim,x:xlim] * 
                                                      mask[y:ylim,x:xlim] )
    
    #Divide the running wavelength by the running weight sum to normalize by
    #how much of the kernel we ended up using
    running_wav_sum[mask] = running_wav_sum[mask]/running_wgt_sum[mask]
    running_wav_sum[np.logical_not(mask)]=wavearray[np.logical_not(mask)]
    
    return running_wav_sum

def convolve_inty(intyarray, badpixarray, kern):
    """Convolves an intensity array by the kernel 'kern'. Ignores pixels with
    values equal to zero so as to preserve user-done masking, chip gaps,
    aperture masking, etc.
    
    Note that the kernel must be square in shape, and furthermore the square
    must be odd in side length, i.e. 1x1, 3x3, 5x5, etc.
    
    Inputs:
    intyarray -> Input intensity array. Masked values ought to be zero.
    badpixarray -> Input bad pixel mask. Ought to be 1 for bad pixels.
    kern -> Kernel for the convolution
    
    Outputs:
    final_inty -> The convolved intensity array
    
    """
    
    #A mask for pixels which are good
    mask = np.logical_and(badpixarray==0,intyarray!=0)
    
    #Kernel size for looping in the most straightforward way
    ksize = (kern.shape[1]-1)/2
    
    #Image dimensions
    xsize = intyarray.shape[1]
    ysize = intyarray.shape[0]
    
    #Empty intensity and count arrays to fill while convolving
    running_int_sum = np.zeros_like(intyarray)
    running_wgt_sum = np.zeros_like(intyarray)
    
    #Loop over the kernel
    for y in range(2*ksize+1):
        for x in range(2*ksize+1):
            ylim = ysize-2*ksize-1+y
            xlim = xsize-2*ksize-1+x
            running_int_sum[ksize:ysize-ksize-1,
                            ksize:xsize-ksize-1] += ( kern[y,x] * 
                                                      intyarray[y:ylim,x:xlim] * 
                                                      mask[y:ylim,x:xlim] )
            running_wgt_sum[ksize:ysize-ksize-1,
                            ksize:xsize-ksize-1] += ( kern[y,x] *
                                                      mask[y:ylim,x:xlim] )
    
    #Divide the running intensity by the running kernel sum to normalize by
    #how much of the kernel we ended up using
    running_int_sum[mask] = running_int_sum[mask]/running_wgt_sum[mask]
    running_int_sum[np.logical_not(mask)]=0
    
    return running_int_sum

def make_final_image(input_image, output_image,
                     desired_fwhm, clobber=False):
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
    desired_fwhm -> Desired FWHM for the resultant image to have.
    
    Optional Inputs:
    clobber -> Overwrite output images if they already exist. Default is False.
    
    """
    
    print "Making final data cube images for image "+input_image
    
    #Open the image
    image = FPImage(input_image)
    
    #Measure the sky background level in the input image
    skyavg, _truesig, skysig = image.skybackground()
    
    #Get various header keywords, crash if necessary
    intygrid = image.inty
    fwhm = image.fwhm
    wave0 = image.wave0
    calf = image.calf
    xcen = image.xcen
    ycen = image.ycen
    if fwhm == None:
        crash("Error! FWHM not measured for image "+input_image+".")
    if wave0 == None or calf == None:
        crash("Error! Wavelength solution does "+
              "not exist for image "+input_image+".")
    if xcen == None or ycen == None:
        crash("Error! Center values not measured "+
              "image "+input_image+".")
    if fwhm>desired_fwhm:
        crash("Error! Desired FWHM too low for image "+
              input_image+".")
    
    #Calculate the necessary FWHM for convolution and make the gaussian kernel
    fwhm_conv = np.sqrt(desired_fwhm**2-fwhm**2)
    sig = fwhm_conv/2.3548+0.0001 #Magic number converts sigma to fwhm
    ksize = np.ceil(4*sig) #Generate the kernel to 4-sigma
    kxgrid, kygrid = np.meshgrid(np.linspace(-ksize,ksize,2*ksize+1),
                                 np.linspace(-ksize,ksize,2*ksize+1))
    kern = np.exp(-(kxgrid**2+kygrid**2)/(2*sig**2)) #Gaussian kernel
    kern = kern/np.sum(kern) #Normalize the kernel
    
    #Convolve the variance array
    vargrid = image.vari
    badpixgrid = image.badp
    
    #Convolve the variances appropriately
    new_vargrid = convolve_variance(vargrid, badpixgrid, kern)
    
    #Add the sky background uncertainty to the uncertainty grid
    vargrid = vargrid+skysig**2
    
    #Create and convolve the wavelength array
    if image.wave is None:
        rgrid = image.rarray(xcen, ycen)
        wavegrid = wave0 / np.sqrt(1+rgrid**2/calf**2)
    else:
        wavegrid = image.wave
    new_wavegrid = convolve_wave(wavegrid, intygrid, badpixgrid, kern)
    
    #Convolve the intensity image
    new_intygrid = convolve_inty(intygrid, badpixgrid, kern)
    new_intygrid[new_intygrid==0] = intygrid[new_intygrid==0]
    image.fwhm = desired_fwhm #Update header FWHM keyword
    
    #Subtract the sky background from the image
    new_intygrid[new_intygrid!=0] -= skyavg
    
    #Create a new fits extension for the wavelength array
    if image.wave is None:
        waveheader = image.openimage[2].header.copy()
        waveheader['EXTVER'] = 4
        wavehdu = ImageHDU(data = new_wavegrid,
                           header = waveheader,
                           name = "WAVE")
        image.openimage[1].header.set('WAVEEXT', 4,
                                      comment='Extension for Wavelength Frame')
        image.openimage.append(wavehdu)
    
    #Put all of the convolved images in the right extensions
    image.inty = new_intygrid
    image.vari = new_vargrid
    image.wave = new_wavegrid
    
    #Write the output file
    image.writeto(output_image, clobber=clobber)
    
    #Close images
    image.close()
    
    return