# This script is not intended to be executed all at once.
# Rather, it records steps needed to run Stan models on HTC

# Resources:
# * How to log in, upload data:
# https://chtc.cs.wisc.edu/uw-research-computing/connecting
# * How to prepare and submit jobs:
# https://chtc.cs.wisc.edu/uw-research-computing/helloworld
# How to run R jobs: https://chtc.cs.wisc.edu/uw-research-computing/r-jobs

# Preliminary: build R packages and Stan
# 
## Locally, transfer build jobs to submit node
scp htc/build-r.sub cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/solanum-aa
scp htc/build-stan.sub cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/solanum-aa

## Login to submit node
ssh cdmuir@ap2002.chtc.wisc.edu

## On submit node
condor_submit -i build-r.sub
tar -xzf R413.tar.gz
export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
R --version # check that R is installed
mkdir packages
export R_LIBS=$PWD/packages
# install packages, quit r
tar -czf packages.tar.gz packages/
# quit interactive job
# R packages are ready
scp cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/solanum-aa/packages.tar.gz htc/ 


condor_submit -i build-stan.sub
tar -xzf cmdstan-2.34.1.tar.gz
cd cmdstan-2.34.1
make build
cd ..
tar -czf cmdstan-2.34.1.tar.gz cmdstan-2.34.1/


# SYNTHETIC DATA

# 1. Compress local directory before transfering
tar --exclude='packages.tar.gz' -czf solanum-aa.tar.gz htc/

# 2. Transfer files before logging in using scp
scp solanum-aa.tar.gz cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir

# 3. Remove tarball locally
rm solanum-aa.tar.gz

# 4. Login to submit node

# 5. Untar on submit node
tar -xzf solanum-aa.tar.gz -C solanum-aa

# 6. Remove tarball from submit node
rm solanum-aa.tar.gz

# 7. Submit jobs
condor_submit solanum-aa/htc/fit_0001.sub # 19164174
condor_submit solanum-aa/htc/fit_0002.sub # 19164175
condor_submit solanum-aa/htc/fit_0003.sub # 19164176

# check status
condor_q 19164174

# 8. Retrieve results
scp cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/fit_sim* objects/

# clean up on submit node
rm fit_*

# ACTUAL DATA

# 1. Transfer files before logging in using scp
scp data/prepared_rh_curves.rds data/stan_rh_curves.rds cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/solanum-aa/htc
scp htc/fit_aa* stan/solanum-aa* cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/solanum-aa/htc

# 2. Login to submit node
ssh cdmuir@ap2002.chtc.wisc.edu

# 3. Submit jobs
condor_submit solanum-aa/htc/fit_aa1.sub # ?
condor_submit solanum-aa/htc/fit_aa2.sub # ?
condor_submit solanum-aa/htc/fit_aa3.sub # ?
condor_submit solanum-aa/htc/fit_aa4.sub # ?
condor_submit solanum-aa/htc/fit_aa5.sub # ?

# check status
# runs with 'automated' scripted model and 4e3 sampling iterations
# first version of automated model with fixed effects of b0, b1, b2
condor_q 25015
condor_q 25016
condor_q 25017
condor_q 25018
condor_q 25019
# runs with 'automated' scripted model and 4e3 sampling iterations
# second version of automated model with random effects of b0, b1, b2


# 4. Retrieve results
scp cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/fit_aa1.rds objects/
scp cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/fit_aa2.rds objects/
scp cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/fit_aa3.rds objects/
scp cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/fit_aa4.rds objects/

# clean up on submit node
rm fit_*
