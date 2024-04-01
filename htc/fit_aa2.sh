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
timeout 10m Rscript fit_aa2.R $1

timeout_exit_status=$?
 
if [ $timeout_exit_status -eq 124 ]; then
    exit 85
fi

exit $timeout_exit_status
