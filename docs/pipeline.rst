Using the Pipeline
==================

Starting Up
-----------

The *saltfppipe* pipeline is designed to be run in the directory above a 'raw/' directory filled with images from the telescope.

You can call it from an interactive environment like python or ipython with::

	>>> from saltfppipe import pipeline
	>>> pipeline()

Or you can call it directly from the command line with::
	
	$ python relative/path/to/saltfppipe/pipeline.py
	
You can also run the pipeline from another location in the directory structure (other than one directory above the 'raw/' image directory::
	
	>>> from saltfppipe import pipeline
	>>> pipeline('relative/path/to/raw/')

Or again from the command line with the path to the raw directory as a calling argument::
	
	$ python relative/path/to/saltfppipe/pipeline.py path/to/raw/)

Creating the 'product/' directory
---------------------------------

The first step of *saltfppipe* is to convert the raw images in 'raw/' to pre-processed images in a 'product/' directory.

If you already have a 'product/' directory, you'll be greeted with the prompt::

	Product directory already exists. Recreate it? (y/n)
	
If you downloaded the 'product/' directory directly from the SALT servers, you should answer **YES** to this question.
The *saltfppipe* and *zSALT* packages create the product directory in a different way than the default SALT pre-processing.

If you have already created a 'product/' directory using *saltfppipe* and don't want to start over, you can answer **NO** to this question.

This step in the reduction process creates a correctly aligned mosaic of SALT RSS's 6 CCD chips, performs crosstalk/gain/bias corrections, and does cosmic ray masking/removal.
It also creates additional image extensions in each image file for keeping track of the noise and bad pixel masks.

.. note::
	I've noticed that occasionally the *zSALT* package gets stuck and hangs forever on my machine during the cosmic ray removal process.
	If you're encountering this problem, locate ``zsalt/imred.py`` and disable multithreading by changing line 103 to read ``multithread=False``.

Image Verification
------------------

Next, the pipeline will display each image in your 'product/' directory and ask you to verify that the image and its header are correct.
This step is not strictly required, but I highly recommend it.
Relatively often, acquisition images of your object may be incorrectly marked as 'ARC' images in their headers.
Or images of Neon/Argon calibration lamps will be marked as 'FLAT' in the headers.

When in the interactive plotting environment, use the 'W' key to flag an image as bad, or the 'S' key to unflag an image.
Use the 'D' and 'A' keys to navigate forward and backward between images.

.. image:: figs/verify/bad.png

This image is an acquisition image, but is marked as 'NGC 2280' in its image headers.
I don't want to use this image, so I've used the 'W' key to flag it as bad. Note the 'Flagged=True' in the plot title.

.. image:: figs/verify/good.png

This image of NGC 2280 is good and its header information appears correct.
Because it's good, I haven't flagged it.

Once the flagging process is completed, any flagged images reappear for a review process.

.. image:: figs/verify/review.png

Here's the acquisition image again. I don't want to use this image, so I press 'D' to delete it.
If its header information had been incorrect, I could have pressed 'H' to input new header information.
If I had flagged it by mistake, I could have pressed 'A' to approve the image.

.. note::
	Pressing 'D' doesn't actually delete files. Rather, it changes a keyword in their headers called 'FPGOOD' to False.
	If you mistakenly 'delete' an image during this process, you can manually edit this header keyword back to True to undo your mistake.
	
Aperture Masking
----------------

SALT FP images are taken with a circular aperture mask in place.
This next step in the pipeline adds pixels outside this circular aperture to the bad pixel mask.
To create the circular aperture, mark a few points near the aperture boundary by clicking on them with the mouse.

.. image:: figs/aperture/blank.png

Here is one of my images of NGC 2280 before I've marked any points along the boundary.
I decide to left click on a point near the bottom right of the image to zoom in on that region.

.. image:: figs/aperture/zoom.png

Here's a zoomed-in view of the area near where I clicked.
I can left click again to mark a point on the boundary, or I can right click to zoom back out.

.. image:: figs/aperture/too_few.png

Now that I've marked three points along the boundary, the code tries to fit a circle to those points.
Note that at the bottom of the image, the fitted circle does not align with the aperture boundary.
Because the points I've marked are all on one side of the image, the circle does not provide a good fit to the aperture.

.. image:: figs/aperture/good.png

After marking a few more points near the aperture all the way around the image, I now have a good fit to the aperture boundary.
Now I right click with my mouse to exit the interactive view.

The pipeline will now add any pixels outside this fitted circle to the bad pixel masks of every image in the 'product/' directory.

Adding DS9 Region Files to the Bad Pixel Mask
---------------------------------------------

Occasionally, you will want to mask a portion of your image that was affected by a satellite trail or some other unexpected phenomenon.
I have added support for this using DS9 region files. If you have an image 'product/image_name.fits' and you wish to mask a geometric portion of it,
create a region file called 'image_name.reg' and place it in the directory *above* the product directory.
In that region file you may have any number of DS9 regions, one per line.

The supported shapes are 'Circle', 'Box', and 'Ellipse'.

The supported file formats are 'DS9/Funtools', 'IRAF PROS', 'CIAO', 'SAOimage', and 'SAOtng'.

When saving the region file, **always use the image coordinate system**; not physical, WCS, or any other.

Measuring Seeing FWHMs
----------------------

Next, the pipeline needs to know the typical full-width at half-max (FWHM) of point sources in your images, also known as the point spread function, or PSF.

.. image:: figs/seeing/blank.png

Another interactive plotting window will appear, this time asking you to press 'Z' to zoom in on a star.
Find a star in your image and place your mouse cursor over it, then press 'Z'.

.. image:: figs/seeing/zoom.png

A zoomed-in plot of the area around this star appears. Again move your cursor over the star, then press 'W' to fit for the star's FWHM.
Alternatively, if you don't like this star, you can press 'Z' to zoom back out.

.. image:: figs/seeing/good.png

The code will now fit a Gaussian profile to that star in each of your images, then display the FWHM as a function of image number.
The mean and standard deviation of stars you've previously fitted will be displayed in red. The newly marked star will be displayed in blue.
To approve the newly marked and fitted star, left click on this plot. To reject it, right click.
The newly marked star in the above plot seems to agree reasonably well with my previously marked stars, so I left click to approve it.

.. image:: figs/seeing/bad.png

This star seems to be systematically blurrier than the several other stars that I've marked, so I right click to reject the fit.

.. image:: figs/seeing/lots.png

Once I've marked several stars and have a fit I'm happy with, I press 'Q' to quit the interactive plotting window.
Alternatively, I can use the 'D' key to delete a star I've previously marked.

After the interactive marking of stars is complete, the mean FWHMs are written to the image heads under the keyword 'FPFWHM'

Ghost Center Fitting
--------------------

Each bright object in a SALT FP image will have a reflection of itself appear at the other side of the image, reflected about a common center.
This center's location is very important later in the data reduction process for finding the wavelength solution.
This next part of the pipeline fits for this center by trying to identify star/ghost pairs in the images.
In dense star fields, this works quite well.
In sparse fields, this routine tends to require some personal intervention.

.. image:: figs/ghosts/good.png

In the dense star field around NGC 2280, the routine has identified dozens of star/ghost pairs and very accurately found the center.
I press 'Q' to approve this fit and move on to the ghost subtraction section.

.. image:: figs/ghosts_manual/bad.png

In this sparse star field around NGC 1325, the routine has failed to find the ghost center, so I press 'W' to flag the center as bad.
The pipeline will not proceed without a good center fit.
When this happens, it's best to try doing things manually.
For that, use the `find_ghost_centers <functions.html#find-ghost-centers>`_ function.

Flat-fielding
-------------

Object Directory
----------------

At this point, the pipeline creates a new directory named after the object you've observed and copies all images associated with that object to the new directory
For the examples I've been using, the pipeline creates a directory called 'NGC2280'.

The purpose of this is to create a 'backup' point before proceeding with the rest of the routines.
Some of these routines are more prone to error, and this preserves a 'pristine' version of the image files in the 'product/' directory.

Image Alignment and Normalization
---------------------------------

Making a Median Image
---------------------

Wavelength Calibration
----------------------

Sky Ring Subtraction
--------------------

Data Cube Creation
------------------

Velocity Map Fitting
--------------------