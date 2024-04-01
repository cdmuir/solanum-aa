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
cd solanum-aa
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
scp checkpoints/*.tar cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/solanum-aa/htc
scp r/functions.R htc/fit_aa* cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/solanum-aa/htc

# 2. Login to submit node
ssh cdmuir@ap2002.chtc.wisc.edu

# 3. Submit jobs
condor_submit -i solanum-aa/htc/fit_aa1.sub # ?
condor_submit solanum-aa/htc/fit_aa1.sub # ?
condor_submit solanum-aa/htc/fit_aa2.sub # ?
condor_submit solanum-aa/htc/fit_aa3.sub # ?
condor_submit solanum-aa/htc/fit_aa4.sub # ?
condor_submit solanum-aa/htc/fit_aa5.sub # ?

# check status
# 4e3 iterations, max_treedepth=12, correct covariance matrix, checkpoint ever 2e2
condor_q 88597 # aa1
condor_q 88598 # aa2
condor_q 88599 # aa3
condor_q 88600 # aa4
condor_q 88601 # aa5
# same as above, but checkpoint every 1e1
condor_q 93809 # aa1
condor_q 93810 # aa2
condor_q 93811 # aa3
condor_q 93812 # aa4
condor_q 93813 # aa5

condor_tail 88597.1
ls /var/lib/condor/spool
vi fit_aa1_40890.log
ls /var/lib/condor/spool/8598/0/cluster88598.proc0.subproc0/chkpt_folder_aa2_0/cp_info
ls /var/lib/condor/spool/8601/0/cluster88601.proc0.subproc0/chkpt_folder_aa5_0/cp_info
cd /home/cdmuir

condor_tail 38668
condor_tail 25015

# 4. Retrieve results
scp cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/draws_aa1_*.rds objects/

scp cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/fit_aa1.rds objects/
scp cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/fit_aa2.rds objects/
scp cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/fit_aa3.rds objects/
scp cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/fit_aa4.rds objects/
scp cdmuir@ap2002.chtc.wisc.edu:/home/cdmuir/fit_aa5.rds objects/

# clean up on submit node
rm fit_*
