#BSUB -n 12
#BSUB -J rstudio
#BSUB -e logs/err.txt
#BSUB -o logs/out.txt
#BSUB -R "rusage[mem=128] span[hosts=1]"
#BSUB -W 08:00
#BSUB -q rna

module load singularity/3.9.2

mkdir -p logs

# path to directory on HPC for persistant storage of R packages
USER_R_LIB=${HOME}/R/catherine/4.2

# What local port to use
LOCAL_PORT=8787

# path to sif file on HPC
SINGULARITY_IMAGE="catherine_bcells2.sif"

# add options for singularity exec
# e.g. "--bind /path/to/some/other/user/directory"
SINGULARITY_EXEC_OPTS=$1

max_n_cores=$(grep processor /proc/cpuinfo | wc -l)
# customize --output path as appropriate (to a directory readable only by the user!)

# Create temporary directory to be populated with directories to bind-mount in the container
# where writable file systems are necessary. Adjust path as appropriate for your computing environment.
workdir=$(python -c 'import tempfile; print(tempfile.mkdtemp())')

mkdir -p -m 700 ${workdir}/run ${workdir}/tmp ${workdir}/var/lib/rstudio-server
cat > ${workdir}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

# Set OMP_NUM_THREADS to prevent OpenBLAS (and any other OpenMP-enhanced
# libraries used by R) from spawning more threads than the number of processors
# allocated to the job.
#
# Set R_LIBS_USER to a path specific to rocker/rstudio to avoid conflicts with
# personal libraries from any R installation in the host environment
cat > ${workdir}/rsession.sh <<END
#!/bin/sh
export OMP_NUM_THREADS=${max_n_cores}
export R_LIBS_USER=$USER_R_LIB
exec /usr/lib/rstudio-server/bin/rsession "\${@}"
END

chmod +x ${workdir}/rsession.sh

export SINGULARITY_BIND="${workdir}/run:/run,${workdir}/tmp:/tmp,${workdir}/database.conf:/etc/rstudio/database.conf,${workdir}/rsession.sh:/etc/rstudio/rsession.sh,${workdir}/var/lib/rstudio-server:/var/lib/rstudio-server"

# Do not suspend idle sessions.
# Alternative to setting session-timeout-minutes=0 in /etc/rstudio/rsession.conf
# https://github.com/rstudio/rstudio/blob/v1.4.1106/src/cpp/server/ServerSessionManager.cpp#L126
export SINGULARITYENV_RSTUDIO_SESSION_TIMEOUT=0

export SINGULARITYENV_USER=$(id -un)
export SINGULARITYENV_PASSWORD=$(openssl rand -base64 15)
# get unused socket per https://unix.stackexchange.com/a/132524
# tiny race condition between the python & singularity commands
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
cat 1>&2 <<END
1. SSH tunnel from your workstation using the following command:
   ssh -N -L $LOCAL_PORT:${HOSTNAME}:${PORT} ${SINGULARITYENV_USER}@LOGIN-HOST
   and point your web browser to http://localhost:$LOCAL_PORT
2. log in to RStudio Server using the following credentials:
   user: ${SINGULARITYENV_USER}
   password: ${SINGULARITYENV_PASSWORD}
When done using RStudio Server, terminate the job by:
1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:
      bkill ${LSB_JOBID}
END

singularity exec $SINGULARITY_EXEC_OPTS --bind /beevol/home/nicholas \
     --cleanenv $SINGULARITY_IMAGE \
    rserver --www-port ${PORT} \
            --server-user ${USER} \
            --auth-none=0 \
            --auth-pam-helper-path=pam-helper \
            --auth-stay-signed-in-days=30 \
            --auth-timeout-minutes=0 \
            --rsession-path=/etc/rstudio/rsession.sh
            
printf 'rserver exited' 1>&2
