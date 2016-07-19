.. saltfppipe documentation master file, created by
   sphinx-quickstart on Tue Jul 12 15:51:10 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*saltfppipe*
============

Welcome to the documentation for *saltfppipe*, a data reduction pipeline for the SALT Fabry-Perot system.

.. note::
	If you make use of this code, please include the following in your acknowledgements section:
	"This research made use of *saltfppipe*, a data reduction package for the SALT Fabry-Perot."
	
	This code is built atop IRAF/Pyraf, PySALT, Astropy, matplotLib, and SciPy. You should also cite/acknowledge these packages as appropriate.

.. note::
	This software was developed primarily to reduce FP data taken in the medium-resolution mode and most of the testing has been with that data.
	It has been used to reduce low-resolution data to some degree of success.
	This code is **not** intended to be used on duel-etalon high-resolution data.

Contents:
---------

.. toctree::
   :maxdepth: 2
   
   install.rst
   pipeline.rst
   fpimage.rst
   functions.rst

This page was last updated by Carl Mitchell on 19 July 2016.