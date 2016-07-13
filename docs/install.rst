Installation
============

Downloading the Package
-----------------------

The latest version of *saltfppipe* is always available from its page on `github <https://github.com/carlmitchell/saltfppipe>`_.
Use the green "Clone or Download" button to either clone the git repository on your local system or download the project as a .zip archive.

Building the Optional *voigtfit* module
---------------------------------------

*saltfppipe* comes with an optional python module called *voigtfit*, which is built from compiled Fortran code.
The module is designed to fit Voigt profiles to the H-alpha line of excited hydrogen and (optionally) the [NII] line 21 Angstroms redward of H-alpha.

If you are using *saltfppipe* to produce H-alpha velocity fields, you'll want to build this module.

If you are using *saltfppipe* only to produce data cubes, the module is unnecessary and you can skip ahead.

In the cloned respository or unzipped directory, execute the following commands:
	* ``$ cd voigtfit``
	* ``$ make``
	* ``$ cd ..``

.. note::
	*voigtfit* requires the *f2py* and *gfortran* packages to be installed

Installing the Package
----------------------

To install the package, run the following command in the cloned repository or unzipped directory:
	* ``$ python setup.py install``

This will install the *saltfppipe* package in your python packages library.

If you built the optional *voigtfit* module, it will be installed automatically.

.. note::
	If your python installation is root-protected, you may need to execute the install command as sudo.

Dependencies
------------

The following is a list of packages on which *saltfppipe* depends. Several of them are either default packages or are readily available from the python package manager "pip":
	* distutils
	* setuptools
	* sys
	* os
	* shutil
	* cStringIO
	* `numpy <http://www.numpy.org/>`_
	* `scipy <https://www.scipy.org/>`_
	* `matplotlib <http://matplotlib.org/>`_
	* `astropy <http://www.astropy.org/>`_
	* `iraf <http://iraf.noao.edu/>`_
	* `pyraf <http://www.stsci.edu/institute/software_hardware/pyraf>`_
	* `zsalt <https://github.com/crawfordsm/zSALT>`_
	* `photutils <https://photutils.readthedocs.io/en/latest/>`_
