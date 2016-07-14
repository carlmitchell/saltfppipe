Useful Functions
================

find_ghost_centers
------------------

In sparse star fields, the default ghost center finding routine occasionally fails.
When this happens, you can use the function "find_ghost_centers" to try to find the center using different input parameters.

.. image:: figs/ghosts_manual/bad.png

In this image of NGC 1325, the default routine has failed to find the reflection center for the ghosts.
I'm pretty sure the center is somewhere around pixel coordinate [792, 502].
We can force it to find this center using the "find_ghost_centers" routine.
This is best done interactively from the python or ipython terminal::

	>>> # The images we want are in the product directory:
	>>> cd product/
	
	>>> # A few useful imports:
	>>> from os import listdir
	>>> from saltfppipe.find_ghost_centers import find_ghost_centers
	
	>>> # I happen to know that the first 3 images and last 1 image
	>>> # of my 'product/' directory are ARC images and not images
	>>> # of NGC 1325, so I don't include them in 'fnlist':
	>>> fnlist = sorted(listdir('.'))[3:-1]
	
	>>> # Now we call the 'find_ghost_centers' function:
	>>> find_ghost_centers(fnlist, tolerance=7,
	>>>                    thresh=3, guess=[792, 502])

.. image:: figs/ghosts_manual/better.png

Here we see that the 'find_ghost_centers' function successfully found several star/ghost pairs and located the center.

The calling keywords are:

	* fnlist -- A list of strings.
	  Each string should be the path to an image file in which you want to detect the ghost center.
	* tolerance -- How close two pixels could be to each other to be considered "close enough", in units of image pixels.
	  Higher values mean the code is more likely to find a center.
	  Extremely high values will likely match objects which aren't truly star/ghost pairs.
	  Default value is 3.
	* thresh -- How bright an object must be in order to be detected, in units of the standard deviation of the sky background.
	  Lower values will detect more objects, increasing the likelihood of finding pairs.
	  Values which are too low will detect spurious objects.
	  Default value is 4.
	* guess -- If you think you know the center, use this keyword.
	  Must be an array-like object of length 2 containing the X and Y coordinates of your guess.
	  No default value.

.. note:: If you're locating pixels in SAOImage DS9, subtract 1 from both coordinates before using it in Python. Python indices begin at zero, while DS9 indices begin at 1.