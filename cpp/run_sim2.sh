# 
# Use the bash shell to interpret this job script 
#$ -S /bin/bash -t 1-100 -l m_mem_free=12G
# 
# Send an e-mail to the address 
# specified in .sge_request when this job ends. 
##$ -m e 
# 
# Only submit this job to nodes 
# that have at least 8GB of RAM free.
##$ -l mem_requested=8G, h_vmem=8G

# Put the hostname, current directory, and start date 
# into variables, then write them to the SGE standard output file.
# Having this information in your output file will help track down
# any errors which might occur during the job's run.
WRKHOST=`/bin/hostname` 
WRKDIR=`/bin/pwd` 
STRTDATE=`/bin/date` 
echo "**** JOB STARTED ON $WRKHOST AT $STRTDATE" 
echo "**** JOB RUNNING IN $WRKDIR" 

# SGE jobs do not run in your login environment, so you'll need to 
# load your environment, or atleast the modules your script needs to run
#source /net/gs/vol1/home/aclark4/.bashrc
# also load any modules as needed
source /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load mpc/latest
module load mpfr/latest
module load gmp/latest
module load gcc/latest
module load gsl/latest

# Script or command(s) to run via SGE

cd ~/go/src/github.com/hyperboliccake/qtl_mapping/cpp/
sleep_time=$((10*$SGE_TASK_ID))
sleep $sleep_time
./sim2_main 1 10000 10000 2 .05 $SGE_TASK_ID
./sim2_main 1 10000 10000 2 .1 $SGE_TASK_ID
./sim2_main 1 10000 10000 2 .01 $SGE_TASK_ID