import numpy as np
from astropy.io.fits import open as openfits
from saltfppipe.wavelib import get_libraries
import matplotlib.pyplot as plt
from saltfppipe.gauss_fit import make_sum_gaussians, nGauss
from scipy.optimize.minpack import curve_fit

class SkyRingPlot:
    def __init__(self, imagearray, subtracted_image_array,
                 xcen, ycen, radii_of_rings,
                 spectrum_x, spectrum_y, spectrum_sub_y,
                 fitted_x, fitted_y,
                 sky_waves, added_waves, extra_waves,
                 imagenumberstring, hasbeensaved, fitsuccess):
        
        self.axis = None
        self.key = None
        self.xcoo = None
        self.ycoo = None
        
        fig = plt.figure()
        
        #Axis 1 - "Before" image plot
        ax1 = fig.add_subplot(221)
        self.plot1 = ax1.imshow(imagearray,
                                cmap = "Greys",
                                aspect = "equal",
                                vmin = np.percentile(imagearray,5),
                                vmax = np.percentile(imagearray,95),
                                origin = "lower",
                                interpolation = "nearest")
        for i in range(len(radii_of_rings)):
            ax1.plot(xcen+radii_of_rings[i]*np.cos(np.linspace(0,2*np.pi,500)),ycen+radii_of_rings[i]*np.sin(np.linspace(0,2*np.pi,500)),color="blue")
        ax1.set_title("Before Image")
        ax1.set_xlim([0,imagearray.shape[1]])
        ax1.set_ylim([0,imagearray.shape[0]])
        plt.tick_params(axis="both",which="both",bottom="off",top="off",right="off",left="off",labelbottom="off",labelleft="off")
        
        #Axis 2 - "Before" spectrum plot
        ax2 = fig.add_subplot(222)
        self.plot2 = ax2.scatter(spectrum_x, spectrum_y,color="blue")
        ax2.plot(fitted_x, fitted_y,color="red")
        for i in range(len(sky_waves)):
            ax2.axvline(sky_waves[i],color="blue")
        for i in range(len(added_waves)):
            ax2.axvline(added_waves[i],color="green")
        for i in range(len(extra_waves)):
            ax2.axvline(extra_waves[i],color='red')
        ax2.set_ylabel("Median Intensity")
        ax2.set_title("Before Spectrum")

        #Axis 3 - "After" image plot
        ax3 = fig.add_subplot(223)
        self.plot3 = ax3.imshow(subtracted_image_array,
                                cmap = "Greys",
                                aspect = "equal",
                                vmin = np.percentile(imagearray,5),
                                vmax = np.percentile(imagearray,95),
                                origin = "lower",
                                interpolation = "nearest")
        #for i in range(len(radii_of_rings)):
        #    ax3.plot(xcen+radii_of_rings[i]*np.cos(np.linspace(0,2*np.pi,500)),ycen+radii_of_rings[i]*np.sin(np.linspace(0,2*np.pi,500)),color="blue")
        ax3.set_xlim([0,imagearray.shape[1]])
        ax3.set_ylim([0,imagearray.shape[0]])
        ax3.set_title("After Image")
        plt.tick_params(axis="both",which="both",bottom="off",top="off",right="off",left="off",labelbottom="off",labelleft="off")

        #Axis 4 - "After" spectrum plot   
        ax4 = fig.add_subplot(224)
        self.plot4 = ax4.scatter(spectrum_x, spectrum_sub_y,color="blue")
        for i in range(len(sky_waves)):
            ax4.axvline(sky_waves[i],color="blue")
        for i in range(len(added_waves)):
            ax4.axvline(added_waves[i],color="green")
        for i in range(len(extra_waves)):
            ax4.axvline(extra_waves[i],color='red')
        ax4.set_xlabel("Wavelength [Ang]")
        ax4.set_ylabel("Median Intensity")
        ax4.set_title("After Spectrum")
        
        if hasbeensaved: savedstring = "already has a saved fit."
        else: savedstring = "does not yet have a saved fit."
        
        if fitsuccess: successstring = "Fit successful!"
        else: successstring = "FIT UNSUCCESSFUL!"
        
        plt.suptitle(successstring+" Image #"+imagenumberstring+" "+savedstring+"\n"+
                     "Fitted lines are blue. Additional known lines are red.\n"+
                     "Press 'a' and 'd' to navigate between images.\n"+
                     "Press 'w' to add a currently-unmarked ring.\n"+
                     "Press 'e' to force-add a ring at the cursor position.\n"+
                     "Press 's' to delete the nearest ring to the cursor.\n"+
                     "Press 'q' to quit and subtract saved profiles from images.\n"+
                     "Press 'r' to save the current ring fit for this image.")
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.68)
                    
        self.cid1 = self.plot1.figure.canvas.mpl_connect('key_press_event',self.keypress)
        self.cid2 = self.plot2.figure.canvas.mpl_connect('key_press_event',self.keypress)
        self.cid3 = self.plot3.figure.canvas.mpl_connect('key_press_event',self.keypress)
        self.cid4 = self.plot4.figure.canvas.mpl_connect('key_press_event',self.keypress)
        
        plt.show()
    
    def keypress(self, event):
        if event.key in ["w","a","s","d","e","q","r"]:
            self.key = event.key
            if event.inaxes in [self.plot1.axes,self.plot2.axes,self.plot3.axes,self.plot4.axes]:
                self.xcoo = event.xdata
                self.ycoo = event.ydata
                if event.inaxes == self.plot1.axes: self.axis = 1
                if event.inaxes == self.plot2.axes: self.axis = 2
                if event.inaxes == self.plot3.axes: self.axis = 3
                if event.inaxes == self.plot4.axes: self.axis = 4
            plt.close()
            return

def sub_sky_rings(fnlist,medfilelist):
    """Fits for night sky rings in a series of images, then subtracts them off.
    If uncertainty images exist, updates them with better uncertainties.
    
    The rings are fitted in an interactive way by the user.
    
    After the ring subtraction, the header keyword "fpdering" is created and
    set to "True."
    
    Note: A 'median.fits' image must exist for the object or this routine will
    not work.
    
    Inputs:
    fnlist -> List of strings, each the path to a fits image.
    medfilelist -> List of paths to median images for each image to subtract off
                   (probably all the same, but generally could be different)
    
    """
    
    #This bit takes care of the 's' to save shortcut in matplotlib.
    oldsavekey = plt.rcParams["keymap.save"]
    plt.rcParams["keymap.save"] = ""
    
    #Open all of the images, median-subtract them, make wavelength arrays
    # make skywavelibs for each, make spectra for each, trim images
    imagelist = []
    wavearraylist = []
    skywavelibs = []
    skywavelibs_extended = []
    spectrum_wave = []
    spectrum_inty = []
    minwavelist = []
    maxwavelist = []
    xcenlist = []
    ycenlist = []
    Flist = []
    wave0list = []
    has_been_saved = np.zeros(len(fnlist),dtype="bool")
    for i in range(len(fnlist)):
        #Open image
        imagelist.append(openfits(fnlist[i]))
        medimage = openfits(medfilelist[i])
        #Get a bunch of relevant header values
        Flist.append(imagelist[i][0].header["fpcalf"])
        wave0list.append(imagelist[i][0].header["fpwave0"])
        xcen = imagelist[i][0].header["fpxcen"]
        ycen = imagelist[i][0].header["fpycen"]
        arad = imagelist[i][0].header["fparad"]
        axcen = imagelist[i][0].header["fpaxcen"]
        aycen = imagelist[i][0].header["fpaycen"]
        #Median subtract it
        imagelist[i][1].data += -medimage[0].data
        #Make wavelength array
        xgrid, ygrid = np.meshgrid(np.arange(imagelist[i][1].data.shape[1]),np.arange(imagelist[i][1].data.shape[0]))
        r2grid = (xgrid-xcen)**2+(ygrid-ycen)**2
        wavearraylist.append(wave0list[i]/np.sqrt(1+r2grid/Flist[i]**2))
        #Get the sky wave library
        minwavelist.append(np.min(wavearraylist[i][r2grid<arad**2]))
        maxwavelist.append(np.max(wavearraylist[i][r2grid<arad**2]))
        skywaves = get_libraries(imagelist[i][0].header["FILTER"])[1]
        skywavelibs.append(skywaves[np.logical_and(skywaves>minwavelist[i],skywaves<maxwavelist[i])])
        #Make an "extended" sky wave library
        skywavelibs_extended.append(skywaves[np.logical_or(np.logical_and(skywaves<minwavelist[i],skywaves>minwavelist[i]-5),np.logical_and(skywaves>maxwavelist[i],skywaves<maxwavelist[i]+5))])
        #Make the spectrum
        spectrum_wave.append(np.linspace(np.min(wavearraylist[i][r2grid<arad**2]),np.max(wavearraylist[i][r2grid<arad**2]),np.int(((np.max(wavearraylist[i][r2grid<arad**2])-np.min(wavearraylist[i][r2grid<arad**2]))/0.25))))
        spectrum_inty.append(np.zeros_like(spectrum_wave[i][1:]))
        for j in range(len(spectrum_wave[i])-1):
            spectrum_inty[i][j] = np.median(imagelist[i][1].data[np.logical_and(np.logical_and(imagelist[i][1].data!=0,wavearraylist[i]>spectrum_wave[i][j]),wavearraylist[i]<spectrum_wave[i][j+1])])
        spectrum_wave[i] = 0.5*(spectrum_wave[i][:-1]+spectrum_wave[i][1:])
        #Trim the images and arrays to the aperture size
        imagelist[i][1].data = imagelist[i][1].data[aycen-arad:aycen+arad,axcen-arad:axcen+arad]
        wavearraylist[i] = wavearraylist[i][aycen-arad:aycen+arad,axcen-arad:axcen+arad]
        xcenlist.append(arad+xcen-axcen)
        ycenlist.append(arad+ycen-aycen)
        #Convert skywavelibs to lists
        skywavelibs[i] = list(skywavelibs[i])
        skywavelibs_extended[i] = list(skywavelibs_extended[i])
        #Close median image
        medimage.close()
    
    #Interactive plotting of ring profiles
    addedwaves = []
    final_fitted_waves = []
    final_fitted_intys = []
    final_fitted_sigs = []
    for i in range(len(fnlist)):
        addedwaves.append([])
        final_fitted_waves.append([])
        final_fitted_intys.append([])
        final_fitted_sigs.append([])
    i = 0
    while True:
        #Determine which wavelengths we want to fit
        waves_to_fit = skywavelibs[i]+addedwaves[i]
        
        #Generate a fitting function from those wavelengths
        func_to_fit = make_sum_gaussians(waves_to_fit)
        
        #Come up with reasonable guesses for the fitting function parameters
        contguess = [0]
        intyguesses = []
        sigguesses = []
        for j in range(len(waves_to_fit)):
            intyguesses.append(np.sum(spectrum_inty[i][np.logical_and(spectrum_wave[i]<waves_to_fit[j]+4,spectrum_wave[i]>waves_to_fit[j]-4)]*(spectrum_wave[i][1]-spectrum_wave[i][0])))
            sigguesses.append(2)
        guess = contguess+intyguesses+sigguesses
        
        #Fit the spectrum
        fitsuccess = True
        try:
            fit = curve_fit(func_to_fit,spectrum_wave[i],spectrum_inty[i],p0=guess)[0]
        except RuntimeError:
            print "Warning: Fit did not converge. Maybe remove some erroneous lines from the fit?"
            fit = guess
            fitsuccess = False
        fitcont = fit[0]
        fitintys = np.array(fit[1:len(waves_to_fit)+1])
        fitsigs = np.array(fit[len(waves_to_fit)+1:2*len(waves_to_fit)+1])
        
        #Make the subtracted plot array
        subarray = imagelist[i][1].data.copy()
        for j in range(len(waves_to_fit)): subarray += -fitintys[j]*nGauss(wavearraylist[i]-waves_to_fit[j],fitsigs[j])
        
        #Figure out which radius each wavelength is at
        radiilist = []
        for j in range(len(waves_to_fit)): radiilist.append(Flist[i]*np.sqrt((wave0list[i]/waves_to_fit[j])**2-1))
        
        #Create the subtracted spectrum
        sub_spec_inty = spectrum_inty[i] - fitcont
        for j in range(len(waves_to_fit)): sub_spec_inty += -fitintys[j]*nGauss(spectrum_wave[i]-waves_to_fit[j],fitsigs[j])
        
        #Create the spectrum fit plot
        fit_X = np.linspace(minwavelist[i],maxwavelist[i],500)
        fit_Y = np.ones_like(fit_X)*fitcont
        for j in range(len(waves_to_fit)): fit_Y += fitintys[j]*nGauss(fit_X-waves_to_fit[j], fitsigs[j])
        
        #Plot the image, spectrum, and fit
        profile_plot = SkyRingPlot(imagelist[i][1].data, subarray,
                                   xcenlist[i], ycenlist[i], radiilist,
                                   spectrum_wave[i], spectrum_inty[i], sub_spec_inty,
                                   fit_X, fit_Y,
                                   skywavelibs[i], addedwaves[i], skywavelibs_extended[i],
                                   repr(i+1)+"/"+repr(len(fnlist)), has_been_saved[i], fitsuccess)
        
        #Shifting images and breakout condition
        if profile_plot.key == "a": i+=-1
        if profile_plot.key == "d": i+=1
        quitloop = False
        if (i==-1 or i==len(fnlist) or profile_plot.key=="q"):
            if (np.sum(np.logical_not(has_been_saved))==0):
                while True:
                    yn = raw_input("Finished fitting ring profiles? (y/n) ")
                    if "n" in yn or "N" in yn:
                        break
                    if "y" in yn or "Y" in yn:
                        quitloop=True
                        break
            else:
                print "Error: Ring fits have not yet been saved for the following images:"
                print (np.arange(len(fnlist))+1)[np.logical_not(has_been_saved)]
        if quitloop: break
        if i == -1: i=0
        if i == len(fnlist): i=len(fnlist)-1
        
        #Delete ring option
        if profile_plot.key == "s":
            nearest_wave = None
            if profile_plot.axis in [1,3]:
                #The click was made in an image plot
                clicked_radius = np.sqrt((profile_plot.xcoo - xcenlist[i])**2 + (profile_plot.ycoo - ycenlist[i])**2)
                nearest_wave = waves_to_fit[np.argmin(np.abs((np.array(radiilist)-clicked_radius)))]
            elif profile_plot.axis in [2,4]:
                #The click was in a spectrum plot
                clicked_wave = profile_plot.xcoo
                nearest_wave = waves_to_fit[np.argmin(np.abs(np.array(waves_to_fit)-clicked_wave))]
            if nearest_wave != None:
                #Remove the nearest wavelength from the sky_wave_lib or added waves, add it to extra waves
                if nearest_wave in skywavelibs[i]:
                    skywavelibs[i].remove(nearest_wave)
                    skywavelibs_extended[i].append(nearest_wave)
                if nearest_wave in addedwaves[i]: addedwaves[i].remove(nearest_wave)
                
        #Add ring option
        if profile_plot.key == "w":
            clicked_wave = None
            if profile_plot.axis in [1,3]:
                #The click was made in an image plot
                clicked_radius = np.sqrt((profile_plot.xcoo - xcenlist[i])**2 + (profile_plot.ycoo - ycenlist[i])**2)
                #Convert radius to wavelength
                clicked_wave = wave0list[i] / (1+clicked_radius**2/Flist[i]**2)**0.5
            elif profile_plot.axis in [2,4]:
                #The click was in a spectrum plot
                clicked_wave = profile_plot.xcoo
            if clicked_wave != None:
                #Find the nearest wavelength in the extended list (if there are any)
                if len(np.abs(np.array(skywavelibs_extended[i])-clicked_wave))>0:
                    wave_to_add = skywavelibs_extended[i][np.argmin(np.abs(np.array(skywavelibs_extended[i])-clicked_wave))]
                    #Add that wave to skywavelibs, remove it from extended
                    skywavelibs[i].append(wave_to_add)
                    skywavelibs_extended[i].remove(wave_to_add)
        
        #Force add ring option
        if profile_plot.key == "e":
            clicked_wave = None
            if profile_plot.axis in [1,3]:
                #The click was made in an image plot
                clicked_radius = np.sqrt((profile_plot.xcoo - xcenlist[i])**2 + (profile_plot.ycoo - ycenlist[i])**2)
                #Convert radius to wavelength
                clicked_wave = wave0list[i] / (1+clicked_radius**2/Flist[i]**2)**0.5
            elif profile_plot.axis in [2,4]:
                #The click was in a spectrum plot
                clicked_wave = profile_plot.xcoo
            if clicked_wave != None:
                #Add this wavelength to the added array
                addedwaves[i].append(clicked_wave)
                
        #Save option
        if profile_plot.key == "r":
            final_fitted_waves[i] = waves_to_fit[:]
            final_fitted_intys[i] = fitintys[:]
            final_fitted_sigs[i] = fitsigs[:]
            has_been_saved[i] = True
        
    #Close all of the images
    medimage.close()
    for i in range(len(fnlist)):
        imagelist[i].close()
    
    #Subtract the ring profiles and update headers
    for i in range(len(fnlist)):
        image = openfits(fnlist[i], mode="update")
        arad = image[0].header["fparad"]
        axcen = image[0].header["fpaxcen"]
        aycen = image[0].header["fpaycen"]
        mask = image[1].data == 0 #For re-correcting chip gaps
        for j in range(len(final_fitted_waves[i])):
            image[1].data[aycen-arad:aycen+arad,axcen-arad:axcen+arad] += -final_fitted_intys[i][j]*nGauss(wavearraylist[i]-final_fitted_waves[i][j], final_fitted_sigs[i][j])
            image[0].header["fpdering"] = "True"
        image[1].data[mask]=0 #For re-correcting chip gaps
        image.close()
    
    #Restore the old keyword shortcut
    plt.rcParams["keymap.save"] = oldsavekey
    
    return