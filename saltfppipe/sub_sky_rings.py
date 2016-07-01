import numpy as np
from astropy.io import fits
from saltfppipe.wavelib import get_libraries
import matplotlib.pyplot as plt
from saltfppipe.gauss_fit import make_sum_gaussians, nGauss
from scipy.optimize.minpack import curve_fit
from saltfppipe.fp_image_class import FPImage


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

        # Axis 1 - "Before" image plot
        ax1 = fig.add_subplot(221)
        self.plot1 = ax1.imshow(imagearray,
                                cmap="Greys",
                                aspect="equal",
                                vmin=np.percentile(imagearray, 5),
                                vmax=np.percentile(imagearray, 95),
                                origin="lower",
                                interpolation="nearest")
        for radius in radii_of_rings:
            ax1.plot(xcen+radius*np.cos(np.linspace(0, 2*np.pi, 500)),
                     ycen+radius*np.sin(np.linspace(0, 2*np.pi, 500)),
                     color="blue")
        ax1.set_title("Before Image")
        ax1.set_xlim([0, imagearray.shape[1]])
        ax1.set_ylim([0, imagearray.shape[0]])
        plt.tick_params(axis="both", which="both",
                        bottom="off", top="off",
                        right="off", left="off",
                        labelbottom="off", labelleft="off")

        # Axis 2 - "Before" spectrum plot
        ax2 = fig.add_subplot(222)
        self.plot2 = ax2.scatter(spectrum_x, spectrum_y, color="blue")
        ax2.plot(fitted_x, fitted_y, color="red")
        for i in range(len(sky_waves)):
            ax2.axvline(sky_waves[i], color="blue")
        for i in range(len(added_waves)):
            ax2.axvline(added_waves[i], color="green")
        for i in range(len(extra_waves)):
            ax2.axvline(extra_waves[i], color='red')
        # Getting minimum and maximum y for plotting spectra
        minyplot = np.min(spectrum_y)
        maxyplot = np.max(spectrum_y)
        yrange = maxyplot - minyplot
        minyplot -= 0.05*yrange
        maxyplot += 0.05*yrange
        ax2.set_ylim(minyplot, maxyplot)

        ax2.set_ylabel("Median Intensity")
        ax2.set_title("Before Spectrum")

        # Axis 3 - "After" image plot
        ax3 = fig.add_subplot(223)
        self.plot3 = ax3.imshow(subtracted_image_array,
                                cmap="Greys",
                                aspect="equal",
                                vmin=np.percentile(imagearray, 5),
                                vmax=np.percentile(imagearray, 95),
                                origin="lower",
                                interpolation="nearest")
        ax3.set_xlim([0, imagearray.shape[1]])
        ax3.set_ylim([0, imagearray.shape[0]])
        ax3.set_title("After Image")
        plt.tick_params(axis="both", which="both",
                        bottom="off", top="off",
                        right="off", left="off",
                        labelbottom="off", labelleft="off")

        # Axis 4 - "After" spectrum plot
        ax4 = fig.add_subplot(224)
        self.plot4 = ax4.scatter(spectrum_x, spectrum_sub_y, color="blue")
        for i in range(len(sky_waves)):
            ax4.axvline(sky_waves[i], color="blue")
        for i in range(len(added_waves)):
            ax4.axvline(added_waves[i], color="green")
        for i in range(len(extra_waves)):
            ax4.axvline(extra_waves[i], color='red')
        # Getting minimum and maximum y for plotting spectra
        minyplot = np.min(spectrum_sub_y)
        maxyplot = np.max(spectrum_sub_y)
        yrange = maxyplot - minyplot
        minyplot -= 0.05*yrange
        maxyplot += 0.05*yrange
        ax4.set_ylim(minyplot, maxyplot)

        ax4.set_xlabel("Wavelength [Ang]")
        ax4.set_ylabel("Median Intensity")
        ax4.set_title("After Spectrum")

        if hasbeensaved:
            savedstring = "already has a saved fit."
        else:
            savedstring = "does not yet have a saved fit."

        if fitsuccess:
            success = "Fit successful!"
        else:
            success = "FIT UNSUCCESSFUL!"

        title = (success+" Image #"+imagenumberstring+" "+savedstring+"\n" +
                 "Fitted lines are blue. Additional known lines are red.\n" +
                 "Press 'a' and 'd' to navigate between images.\n" +
                 "Press 'w' to add a currently-unmarked ring.\n" +
                 "Press 'e' to force-add a ring at the cursor position.\n" +
                 "Press 's' to delete the nearest ring to the cursor.\n" +
                 "Press 'q' to quit and subtract saved profiles.\n" +
                 "Press 'r' to save the current ring fit for this image.\n" +
                 "Press 'x' to switch to manual mode.")

        plt.suptitle(title)

        plt.tight_layout()
        plt.subplots_adjust(top=0.65)

        self.cid1 = self.plot1.figure.canvas.mpl_connect('key_press_event',
                                                         self.keypress)
        self.cid2 = self.plot2.figure.canvas.mpl_connect('key_press_event',
                                                         self.keypress)
        self.cid3 = self.plot3.figure.canvas.mpl_connect('key_press_event',
                                                         self.keypress)
        self.cid4 = self.plot4.figure.canvas.mpl_connect('key_press_event',
                                                         self.keypress)

        plt.show()

    def keypress(self, event):
        if event.key in "wasdeqrx":
            self.key = event.key
            if event.inaxes in [self.plot1.axes, self.plot2.axes,
                                self.plot3.axes, self.plot4.axes]:
                self.xcoo = event.xdata
                self.ycoo = event.ydata
                if event.inaxes == self.plot1.axes:
                    self.axis = 1
                if event.inaxes == self.plot2.axes:
                    self.axis = 2
                if event.inaxes == self.plot3.axes:
                    self.axis = 3
                if event.inaxes == self.plot4.axes:
                    self.axis = 4
            plt.close()
            return


class ManualPlot:
    def __init__(self,
                 imagearray, subtracted_image_array,
                 xcen, ycen,
                 spectrum_x, spectrum_y, spectrum_sub_y,
                 fitted_x, fitted_y,
                 hasbeensaved, imagenumberstring,
                 radii_of_rings, wavelengths_of_lines,
                 colors):

        self.axis = None
        self.key = None
        self.xcoo = None
        self.ycoo = None

        fig = plt.figure()

        # Axis 1 - "Before" image plot
        ax1 = fig.add_subplot(221)
        self.plot1 = ax1.imshow(imagearray,
                                cmap="Greys",
                                aspect="equal",
                                vmin=np.percentile(imagearray, 5),
                                vmax=np.percentile(imagearray, 95),
                                origin="lower",
                                interpolation="nearest")
        for (r, c) in zip(radii_of_rings, colors):
            if not np.isnan(r):
                ax1.plot(xcen+r*np.cos(np.linspace(0, 2*np.pi, 500)),
                         ycen+r*np.sin(np.linspace(0, 2*np.pi, 500)),
                         color=c)
        ax1.set_title("Before Image")
        ax1.set_xlim([0, imagearray.shape[1]])
        ax1.set_ylim([0, imagearray.shape[0]])
        plt.tick_params(axis="both", which="both",
                        bottom="off", top="off",
                        right="off", left="off",
                        labelbottom="off", labelleft="off")

        # Axis 2 - "Before" spectrum plot
        ax2 = fig.add_subplot(222)
        self.plot2 = ax2.scatter(spectrum_x, spectrum_y, color="blue")
        ax2.plot(fitted_x, fitted_y, color="red")
        for i in range(len(wavelengths_of_lines)):
            ax2.axvline(wavelengths_of_lines[i], color=colors[i])

        # Getting minimum and maximum y for plotting spectra
        minyplot = np.min(spectrum_y)
        maxyplot = np.max(spectrum_y)
        yrange = maxyplot - minyplot
        minyplot -= 0.05*yrange
        maxyplot += 0.05*yrange
        ax2.set_ylim(minyplot, maxyplot)

        ax2.set_ylabel("Median Intensity")
        ax2.set_title("Before Spectrum")

        # Axis 3 - "After" image plot
        ax3 = fig.add_subplot(223)
        self.plot3 = ax3.imshow(subtracted_image_array,
                                cmap="Greys",
                                aspect="equal",
                                vmin=np.percentile(imagearray, 5),
                                vmax=np.percentile(imagearray, 95),
                                origin="lower",
                                interpolation="nearest")
        ax3.set_xlim([0, imagearray.shape[1]])
        ax3.set_ylim([0, imagearray.shape[0]])
        ax3.set_title("After Image")
        plt.tick_params(axis="both", which="both",
                        bottom="off", top="off",
                        right="off", left="off",
                        labelbottom="off", labelleft="off")

        # Axis 4 - "After" spectrum plot
        ax4 = fig.add_subplot(224)
        self.plot4 = ax4.scatter(spectrum_x, spectrum_sub_y, color="blue")
        for i in range(len(wavelengths_of_lines)):
            ax4.axvline(wavelengths_of_lines[i], color=colors[i])

        # Getting minimum and maximum y for plotting spectra
        minyplot = np.min(spectrum_sub_y)
        maxyplot = np.max(spectrum_sub_y)
        yrange = maxyplot - minyplot
        minyplot -= 0.05*yrange
        maxyplot += 0.05*yrange
        ax4.set_ylim(minyplot, maxyplot)

        ax4.set_xlabel("Wavelength [Ang]")
        ax4.set_ylabel("Median Intensity")
        ax4.set_title("After Spectrum")

        if hasbeensaved:
            savedstring = "already has a saved fit."
        else:
            savedstring = "does not yet have a saved fit."

        title = ("Image #"+imagenumberstring+" "+savedstring+"\n" +
                 "Included lines are blue. Additional known lines are red.\n" +
                 "Added lines are green. Selected line is purple.\n" +
                 "Press 'z' to cycle selected line, 'r' to save fit.\n" +
                 "Press 'c/v' to make refinements more/less fine.\n" +
                 "Press 't/g' to increase/decrease continuum strength.\n" +
                 "Press 'w/s' to increase/decrease selected line strength.\n" +
                 "Press 'q/a' to increase/decrease selected line width.\n" +
                 "Press 'e/d' to add/remove known lines, " +
                 "or 'f' to force-add.\n" +
                 "Press 'x' to switch to automatic mode (and move on).")

        plt.suptitle(title)

        plt.tight_layout()
        plt.subplots_adjust(top=0.62)

        self.cid1 = self.plot1.figure.canvas.mpl_connect('key_press_event',
                                                         self.keypress)
        self.cid2 = self.plot2.figure.canvas.mpl_connect('key_press_event',
                                                         self.keypress)
        self.cid3 = self.plot3.figure.canvas.mpl_connect('key_press_event',
                                                         self.keypress)
        self.cid4 = self.plot4.figure.canvas.mpl_connect('key_press_event',
                                                         self.keypress)
        plt.show()

    def keypress(self, event):
        if event.key in "zrcvwsqaefdxtg":
            self.key = event.key
            if event.inaxes in [self.plot1.axes, self.plot2.axes,
                                self.plot3.axes, self.plot4.axes]:
                self.xcoo = event.xdata
                self.ycoo = event.ydata
                if event.inaxes == self.plot1.axes:
                    self.axis = 1
                if event.inaxes == self.plot2.axes:
                    self.axis = 2
                if event.inaxes == self.plot3.axes:
                    self.axis = 3
                if event.inaxes == self.plot4.axes:
                    self.axis = 4
            plt.close()
            return


def sub_sky_rings(fnlist, medfilelist):
    """Fits for night sky rings in a series of images, then subtracts them off.
    If uncertainty images exist, updates them with better uncertainties.

    The rings are fitted in an interactive way by the user.

    After the ring subtraction, the header keyword "fpdering" is created and
    set to "True."

    Note: A 'median.fits' image must exist for the object or this routine will
    not work.

    Inputs:
    fnlist -> List of strings, each the path to a fits image.
    medfilelist -> List of paths to median images for each image to subtract
                   off (probably all the same, but could be different)

    """

    # This bit takes care of the 's' to save shortcut in matplotlib.
    oldsavekey = plt.rcParams["keymap.save"]
    plt.rcParams["keymap.save"] = ""
    oldfullscreenkey = plt.rcParams['keymap.fullscreen']
    plt.rcParams['keymap.fullscreen'] = ""

    # Open all of the images, median-subtract them, make wavelength arrays
    # make skywavelibs for each, make spectra for each, trim images
    images = [FPImage(fn) for fn in fnlist]
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
    has_been_saved = np.zeros(len(fnlist), dtype="bool")
    for i in range(len(fnlist)):
        # Get a bunch of relevant header values
        medimage = fits.open(medfilelist[i])
        Flist.append(images[i].calf)
        wave0list.append(images[i].wave0)
        xcen = images[i].xcen
        ycen = images[i].ycen
        arad = images[i].arad
        axcen = images[i].axcen
        aycen = images[i].aycen
        # Median subtract it
        images[i].inty -= medimage[0].data
        images[i].inty -= np.median(images[i].inty[images[i].badp != 1])
        # Make wavelength array
        r2grid = images[i].rarray(xcen, ycen)**2
        wavearraylist.append(wave0list[i]/np.sqrt(1+r2grid/Flist[i]**2))
        # Get the sky wave library
        minwavelist.append(np.min(wavearraylist[i][r2grid < arad**2]))
        maxwavelist.append(np.max(wavearraylist[i][r2grid < arad**2]))
        skywaves = get_libraries(images[i].filter)[1]
        skywavelibs.append(skywaves[np.logical_and(skywaves > minwavelist[i],
                                                   skywaves < maxwavelist[i])])
        # Make an "extended" sky wave library
        lowermask = np.logical_and(skywaves < minwavelist[i],
                                   skywaves > minwavelist[i]-5)
        uppermask = np.logical_and(skywaves > maxwavelist[i],
                                   skywaves < maxwavelist[i]+5)
        skywavelibs_extended.append(skywaves[np.logical_or(uppermask,
                                                           lowermask)])
        # Make the spectrum
        minwav = np.min(wavearraylist[i][r2grid < arad**2])
        maxwav = np.max(wavearraylist[i][r2grid < arad**2])
        wavstep = 0.25  # Quarter-angstrom spacing
        spectrum_wave.append(np.linspace(minwav, maxwav,
                                         np.int(((maxwav-minwav)/wavstep))))
        spectrum_inty.append(np.zeros_like(spectrum_wave[i][1:]))
        for j in range(len(spectrum_wave[i])-1):
            uppermask = wavearraylist[i] < spectrum_wave[i][j+1]
            lowermask = wavearraylist[i] > spectrum_wave[i][j]
            mask = np.logical_and(uppermask, lowermask)
            goodmask = np.logical_and(mask, images[i].badp == 0)
            spectrum_inty[i][j] = np.median(images[i].inty[goodmask])
        spectrum_wave[i] = 0.5*(spectrum_wave[i][:-1]+spectrum_wave[i][1:])
        # Trim the images and arrays to the aperture size
        images[i].inty = images[i].inty[aycen-arad:aycen+arad,
                                        axcen-arad:axcen+arad]
        wavearraylist[i] = wavearraylist[i][aycen-arad:aycen+arad,
                                            axcen-arad:axcen+arad]
        xcenlist.append(arad+xcen-axcen)
        ycenlist.append(arad+ycen-aycen)
        # Convert skywavelibs to lists
        skywavelibs[i] = list(skywavelibs[i])
        skywavelibs_extended[i] = list(skywavelibs_extended[i])
        # Close median image
        medimage.close()
        # Try to fix the NaNs
        mask = np.logical_not(np.isnan(spectrum_inty[i]))
        spectrum_wave[i] = spectrum_wave[i][mask]
        spectrum_inty[i] = spectrum_inty[i][mask]

    # Interactive plotting of ring profiles
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
        # Determine which wavelengths we want to fit
        waves_to_fit = skywavelibs[i]+addedwaves[i]

        # Generate a fitting function from those wavelengths
        func_to_fit = make_sum_gaussians(waves_to_fit)

        # Come up with reasonable guesses for the fitting function parameters
        contguess = [0]
        intyguesses = []
        sigguesses = []
        for j in range(len(waves_to_fit)):
            uppermask = spectrum_wave[i] < waves_to_fit[j]+4
            lowermask = spectrum_wave[i] > waves_to_fit[j]-4
            mask = np.logical_and(uppermask, lowermask)
            # Guess at line intensity via integral
            intyguesses.append(np.sum(spectrum_inty[i][mask] *
                                      (spectrum_wave[i][1] -
                                       spectrum_wave[i][0])))
            # Guess sigma is 2 angstroms
            sigguesses.append(2)
        guess = contguess+intyguesses+sigguesses

        # Fit the spectrum
        fitsuccess = True
        try:
            fit = curve_fit(func_to_fit, spectrum_wave[i],
                            spectrum_inty[i], p0=guess)[0]
        except RuntimeError:
            print ("Warning: Fit did not converge. " +
                   "Maybe remove some erroneous lines from the fit?")
            fit = guess
            fitsuccess = False
        fitcont = fit[0]
        fitintys = np.array(fit[1:len(waves_to_fit)+1])
        fitsigs = np.array(fit[len(waves_to_fit)+1:2*len(waves_to_fit)+1])
        # Flip signs of negative sigma fits (sign degeneracy)
        fitintys[fitsigs < 0] = -1*fitintys[fitsigs < 0]
        fitsigs[fitsigs < 0] = -1*fitsigs[fitsigs < 0]

        # Make the subtracted plot array
        subarray = images[i].inty.copy()
        for j in range(len(waves_to_fit)):
            subarray -= fitintys[j]*nGauss(wavearraylist[i]-waves_to_fit[j],
                                           fitsigs[j])

        # Figure out which radius each wavelength is at
        radiilist = []
        for j in range(len(waves_to_fit)):
            radiilist.append(Flist[i] *
                             np.sqrt((wave0list[i]/waves_to_fit[j])**2-1))

        # Create the subtracted spectrum
        sub_spec_inty = spectrum_inty[i] - fitcont
        for j in range(len(waves_to_fit)):
            sub_spec_inty -= (fitintys[j] *
                              nGauss(spectrum_wave[i]-waves_to_fit[j],
                                     fitsigs[j]))

        # Create the spectrum fit plot
        fit_X = np.linspace(minwavelist[i], maxwavelist[i], 500)
        fit_Y = np.ones_like(fit_X)*fitcont
        for j in range(len(waves_to_fit)):
            fit_Y += fitintys[j]*nGauss(fit_X-waves_to_fit[j], fitsigs[j])

        # Plot the image, spectrum, and fit
        profile_plot = SkyRingPlot(images[i].inty, subarray,
                                   xcenlist[i], ycenlist[i], radiilist,
                                   spectrum_wave[i], spectrum_inty[i],
                                   sub_spec_inty, fit_X, fit_Y,
                                   skywavelibs[i], addedwaves[i],
                                   skywavelibs_extended[i],
                                   repr(i+1)+"/"+repr(len(fnlist)),
                                   has_been_saved[i], fitsuccess)

        # Shifting images and breakout condition
        if profile_plot.key == "a":
            i += -1
        if profile_plot.key == "d":
            i += 1
        quitloop = False
        if (i == -1 or i == len(fnlist) or profile_plot.key == "q"):
            if (np.sum(np.logical_not(has_been_saved)) == 0):
#                 while True:
#                     yn = raw_input("Finished fitting ring profiles? (y/n) ")
#                     if "n" in yn or "N" in yn:
#                         break
#                     if "y" in yn or "Y" in yn:
#                         quitloop = True
#                         break
                quitloop = True
            else:
                print ("Error: Ring fits have not yet been " +
                       "saved for the following images:")
                imagenumbers = np.arange(len(fnlist))+1
                print imagenumbers[np.logical_not(has_been_saved)]
        if quitloop:
            break
        if i == -1:
            i = 0
        if i == len(fnlist):
            i = len(fnlist)-1

        # Delete ring option
        if profile_plot.key == "s":
            nearest_wave = None
            if profile_plot.axis in [1, 3] and len(waves_to_fit) != 0:
                # The keypress was made in an image plot
                clicked_radius = np.sqrt((profile_plot.xcoo - xcenlist[i])**2 +
                                         (profile_plot.ycoo - ycenlist[i])**2)
                near_index = np.argmin(np.abs((np.array(radiilist) -
                                               clicked_radius)))
                nearest_wave = waves_to_fit[near_index]
            elif profile_plot.axis in [2, 4] and len(waves_to_fit) != 0:
                # The keypress was in a spectrum plot
                clicked_wave = profile_plot.xcoo
                near_index = np.argmin(np.abs(np.array(waves_to_fit) -
                                              clicked_wave))
                nearest_wave = waves_to_fit[near_index]
            if nearest_wave is not None:
                # Remove the nearest wavelength from the sky_wave_lib or
                # added waves, add it to extra waves
                if nearest_wave in skywavelibs[i]:
                    skywavelibs[i].remove(nearest_wave)
                    skywavelibs_extended[i].append(nearest_wave)
                if nearest_wave in addedwaves[i]:
                    addedwaves[i].remove(nearest_wave)

        # Add ring option
        if profile_plot.key == "w":
            clicked_wave = None
            if profile_plot.axis in [1, 3]:
                # The keypress was made in an image plot
                clicked_radius = np.sqrt((profile_plot.xcoo - xcenlist[i])**2 +
                                         (profile_plot.ycoo - ycenlist[i])**2)
                # Convert radius to wavelength
                clicked_wave = (wave0list[i] /
                                (1+clicked_radius**2/Flist[i]**2)**0.5)
            elif profile_plot.axis in [2, 4]:
                # The keypress was in a spectrum plot
                clicked_wave = profile_plot.xcoo
            if clicked_wave is not None:
                # Find the nearest wavelength in the extended list
                if len(np.array(skywavelibs_extended[i])-clicked_wave) > 0:
                    near = np.argmin(np.abs(np.array(skywavelibs_extended[i]) -
                                            clicked_wave))
                    wave_to_add = skywavelibs_extended[i][near]
                    # Add that wave to skywavelibs, remove it from extended
                    skywavelibs[i].append(wave_to_add)
                    skywavelibs_extended[i].remove(wave_to_add)

        # Force add ring option
        if profile_plot.key == "e":
            clicked_wave = None
            if profile_plot.axis in [1, 3]:
                # The click was made in an image plot
                clicked_radius = np.sqrt((profile_plot.xcoo - xcenlist[i])**2 +
                                         (profile_plot.ycoo - ycenlist[i])**2)
                # Convert radius to wavelength
                clicked_wave = (wave0list[i] /
                                (1+clicked_radius**2/Flist[i]**2)**0.5)
            elif profile_plot.axis in [2, 4]:
                # The click was in a spectrum plot
                clicked_wave = profile_plot.xcoo
            if clicked_wave is not None:
                # Add this wavelength to the added array
                addedwaves[i].append(clicked_wave)

        # Save option
        if profile_plot.key == "r":
            final_fitted_waves[i] = waves_to_fit[:]
            final_fitted_intys[i] = fitintys[:]
            final_fitted_sigs[i] = fitsigs[:]
            has_been_saved[i] = True

        # Manual Mode Option
        if profile_plot.key == "x":

            # Initialize the manual mode lists
            man_waves = list(waves_to_fit)
            man_intys = list(fitintys)
            man_sigs = list(fitsigs)
            man_cont = fitcont
            man_real = []
            for wave in man_waves:
                if wave in skywavelibs[i] or wave in skywavelibs_extended[i]:
                    man_real.append(True)
                else:
                    man_real.append(False)

            # Initialize refinement level and selected line
            selected_line = 0
            refine = 1
            inty_refine = 0.25*(np.max(spectrum_inty[i]) -
                                np.min(spectrum_inty[i]))
            sig_refine = 0.5
            cont_refine = inty_refine

            # Loop until completion
            while True:

                # Make the subtracted plot array
                subarray = images[i].inty.copy()
                for j in range(len(man_waves)):
                    subarray -= man_intys[j]*nGauss((wavearraylist[i] -
                                                     man_waves[j]),
                                                    man_sigs[j])

                # Figure out which radius each wavelength is at
                radiilist = []
                for j in range(len(man_waves)):
                    radiilist.append(Flist[i] *
                                     np.sqrt((wave0list[i]/man_waves[j])**2-1))
                for skywav in skywavelibs_extended[i]:
                    radiilist.append(Flist[i] *
                                     np.sqrt((wave0list[i] / skywav)**2-1))

                # Create the subtracted spectrum
                sub_spec_inty = spectrum_inty[i] - man_cont
                for j in range(len(man_waves)):
                    sub_spec_inty -= (man_intys[j] *
                                      nGauss(spectrum_wave[i]-man_waves[j],
                                             man_sigs[j]))

                # Create the spectrum fit plot
                fit_X = np.linspace(minwavelist[i], maxwavelist[i], 500)
                fit_Y = np.ones_like(fit_X)*man_cont
                for j in range(len(man_waves)):
                    fit_Y += man_intys[j]*nGauss(fit_X-man_waves[j],
                                                 man_sigs[j])

                # What colors are we using?
                colors = []
                for k in range(len(man_waves)):
                    if selected_line == k:
                        colors.append('purple')
                    elif not man_real[k]:
                        colors.append('green')
                    elif man_real[k]:
                        colors.append('blue')
                for _w in skywavelibs_extended[i]:
                    colors.append('red')

                # Plot the spectrum and fit, get user input
                man_plot = ManualPlot(images[i].inty, subarray,
                                      xcenlist[i], ycenlist[i],
                                      spectrum_wave[i], spectrum_inty[i],
                                      sub_spec_inty,
                                      fit_X, fit_Y,
                                      has_been_saved[i], (repr(i+1)+"/" +
                                                          repr(len(fnlist))),
                                      radiilist,
                                      man_waves+skywavelibs_extended[i],
                                      colors)

                # Cycle selected line
                if man_plot.key == "z":
                    selected_line += 1
                    if selected_line >= len(man_waves):
                        selected_line = 0
                    refine = 1

                # Save fit
                if man_plot.key == "r":
                    final_fitted_waves[i] = man_waves[:]
                    final_fitted_intys[i] = man_intys[:]
                    final_fitted_sigs[i] = man_sigs[:]
                    has_been_saved[i] = True
                    break

                # Refinement levels
                if man_plot.key == "c":
                    refine *= 0.5
                if man_plot.key == "v":
                    refine *= 2

                # Increase/decrease inty
                if len(man_waves) != 0:
                    if man_plot.key == "w":
                        man_intys[selected_line] += refine*inty_refine
                    if man_plot.key == "s":
                        man_intys[selected_line] -= refine*inty_refine

                # Increase/decrease line width
                if len(man_waves) != 0:
                    if man_plot.key == "q":
                        man_sigs[selected_line] += refine*sig_refine
                    if man_plot.key == "a":
                        man_sigs[selected_line] -= refine*sig_refine

                # Increase/decrease continuum
                if man_plot.key == "t":
                    man_cont += refine*cont_refine
                if man_plot.key == "g":
                    man_cont -= refine*cont_refine

                # Add the known line nearest the keypress
                if man_plot.key == "e":
                    clicked_wave = None
                    if man_plot.axis in [1, 3]:
                        # The keypress was made in an image plot
                        clicked_radius = np.sqrt((man_plot.xcoo -
                                                  xcenlist[i])**2 +
                                                 (man_plot.ycoo -
                                                  ycenlist[i])**2)
                        # Convert radius to wavelength
                        clicked_wave = (wave0list[i] /
                                        (1+clicked_radius**2/Flist[i]**2)**0.5)
                    elif man_plot.axis in [2, 4]:
                        # The keypress was in a spectrum plot
                        clicked_wave = man_plot.xcoo
                    if clicked_wave is not None:
                        # Find the nearest wavelength in the extended list
                        if len(np.array(skywavelibs_extended[i]) -
                               clicked_wave) > 0:
                            sky = np.array(skywavelibs_extended[i])
                            near = np.argmin(np.abs(sky-clicked_wave))
                            wave_to_add = skywavelibs_extended[i][near]
                            # Add that wave to manwaves,
                            # remove from extended list
                            man_waves.append(wave_to_add)
                            man_intys.append(0)
                            man_sigs.append(1)
                            man_real.append(True)
                            skywavelibs_extended[i].remove(wave_to_add)
                            selected_line = len(man_waves)-1

                # Remove line nearest the keypress
                if man_plot.key == "d":
                    nearest_wave = None
                    if man_plot.axis in [1, 3] and len(man_waves) != 0:
                        # The keypress was made in an image plot
                        clicked_radius = np.sqrt((man_plot.xcoo -
                                                  xcenlist[i])**2 +
                                                 (man_plot.ycoo -
                                                  ycenlist[i])**2)
                        near_index = np.argmin(np.abs((np.array(radiilist) -
                                                       clicked_radius)))
                        nearest_wave = man_waves[near_index]
                    elif man_plot.axis in [2, 4] and len(man_waves) != 0:
                        # The keypress was in a spectrum plot
                        clicked_wave = man_plot.xcoo
                        near_index = np.argmin(np.abs(np.array(man_waves) -
                                                      clicked_wave))
                        nearest_wave = man_waves[near_index]
                    if nearest_wave is not None:
                        # Remove the nearest wavelength
                        # re-add to extended list
                        # if it was real
                        wave_index = np.where(np.array(man_waves) ==
                                              nearest_wave)[0][0]
                        man_waves.pop(wave_index)
                        man_intys.pop(wave_index)
                        man_sigs.pop(wave_index)
                        if man_real[wave_index]:
                            skywavelibs_extended[i].append(nearest_wave)
                        man_real.pop(wave_index)
                        if wave_index == selected_line:
                            selected_line = 0

                # Force-add
                if man_plot.key == "f":
                    clicked_wave = None
                    if profile_plot.axis in [1, 3]:
                        # The click was made in an image plot
                        clicked_radius = np.sqrt((man_plot.xcoo -
                                                  xcenlist[i])**2 +
                                                 (man_plot.ycoo -
                                                  ycenlist[i])**2)
                        # Convert radius to wavelength
                        clicked_wave = (wave0list[i] /
                                        (1+clicked_radius**2 /
                                         Flist[i]**2)**0.5)
                    elif profile_plot.axis in [2, 4]:
                        # The click was in a spectrum plot
                        clicked_wave = man_plot.xcoo
                    if clicked_wave is not None:
                        # Add this wavelength to the added array
                        man_waves.append(clicked_wave)
                        man_intys.append(0)
                        man_sigs.append(1)
                        man_real.append(False)
                        selected_line = len(man_waves)-1

                # Switch back to automatic
                if man_plot.key == "x":
                    break

    # Close all of the images
    for image in images:
        image.close()

    # Subtract the ring profiles and update headers
    for i in range(len(fnlist)):
        image = FPImage(fnlist[i], update=True)
        arad = image.arad
        axcen = image.axcen
        aycen = image.aycen
        mask = image.inty == 0  # For re-correcting chip gaps
        for j in range(len(final_fitted_waves[i])):
            image.header['subwave'+repr(j)] = final_fitted_waves[i][j]
            image.header['subinty'+repr(j)] = final_fitted_intys[i][j]
            image.header['subsig'+repr(j)] = final_fitted_sigs[i][j]
            ringimg = (final_fitted_intys[i][j] *
                       nGauss(wavearraylist[i]-final_fitted_waves[i][j],
                              final_fitted_sigs[i][j]))
            image.inty[aycen-arad:aycen+arad, axcen-arad:axcen+arad] -= ringimg
            image.ringtog = "True"
        image.inty[mask] = 0  # For re-correcting chip gaps
        image.close()

    # Restore the old keyword shortcut
    plt.rcParams["keymap.save"] = oldsavekey
    plt.rcParams['keymap.fullscreen'] = oldfullscreenkey

    return
