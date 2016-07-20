from numpy.distutils.core import setup
from os.path import isfile, join


# Janky way to install the voigtfit module
if isfile(join('voigtfit', 'voigtfit.so')):
    df = [('', [join('voigtfit', 'voigtfit.so')])]
else:
    df = []

setup(name='saltfppipe',
      version='0.0.1',
      description='SALT Fabry-Perot Data Reduction Software',
      author='Carl J. Mitchell',
      author_email='cmitchell@physics.rutgers.edu',
      packages=['saltfppipe'],
      data_files=df
      )
