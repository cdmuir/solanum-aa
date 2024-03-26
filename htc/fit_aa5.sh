#!/bin/bash

# untar R installation
tar -xzf R413.tar.gz

# untar Stan installation
tar -xzf cmdstan-2.34.1.tar.gz

# untar R Packages
tar -xzf packages.tar.gz

# make sure the script will use your R installation, 
# and the working directory as its home location
export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
export R_LIBS=$PWD/packages

# run your script
Rscript fit_aa5.R $1
