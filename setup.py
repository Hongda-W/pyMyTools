from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='pyMyTools',
      version='0.0.1',
      description='Small Python tools',
      long_description=long_description,
      long_description_content_type="text/markdown",
      packages=['pymytools'],
      author='Hongda Wang',
      author_email='hongda.wang@colorado.edu',
      url='https://github.com/Hongda-W',
      include_package_data=True,
      install_requires= [
          "numpy",
          "matplotlib",
          "Basemap",
          "obspy",
          "configparser",
          "pycpt",
          "netCDF4",
          ], 
      classifiers=[
          'Development Status :: 1 - Planning',
          "License :: OSI Approved :: MIT License",
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3',
          ],
     )
