from distutils.core import setup
from setuptools import find_packages
from os.path import isfile, join
from shutil import copyfile

if isfile(join('voigtfit', 'voigtfit.so')):
    copyfile(join('voigtfit', 'voigtfit.so'),
             join('saltfppipe', 'voigtfit.so'))

setup(name='saltfppipe',
      version='0.0.1',
      description='SALT Fabry-Perot Data Reduction Software',
      author='Carl J. Mitchell',
      author_email='cmitchell@physics.rutgers.edu',
      packages=find_packages()
      )
