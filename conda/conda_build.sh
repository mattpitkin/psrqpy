#!/bin/bash

# This script builds the conda package for psrqpy from the existing PyPI version
# of the package. It is adapted from the script decribed here
# https://medium.com/@giswqs/building-a-conda-package-and-uploading-it-to-anaconda-cloud-6a3abd1c5c52
# and found here https://gist.github.com/giswqs/4eb62fb08658c8a200c4e18bb5e6270c

# The package name from PyPi
pkg='psrqpy'

# the Python versions to build
array=( 2.7 3.5 3.6 3.7 )
echo "Building conda package ..."

# create meta.yaml package information file
conda skeleton pypi $pkg

# building conda packages
for i in "${array[@]}"
do
	conda-build --python $i $pkg
done

# convert package to other platforms
#platforms=( linux-64 )
find $CONDA_PREFIX/conda-bld/linux-64/ -name ${pkg}*.tar.bz2 | while read file
do
    echo $file
    conda convert --platform all $file  -o $CONDA_PREFIX/conda-bld/
    #for platform in "${platforms[@]}"
    #do
    #   conda convert --platform $platform $file  -o $HOME/conda-bld/
    #done
done

# upload packages to conda
find $CONDA_PREFIX/conda-bld/ -name ${pkg}*.tar.bz2 | while read file
do
    echo $file
    anaconda upload $file
done
echo "Building conda package done!"

