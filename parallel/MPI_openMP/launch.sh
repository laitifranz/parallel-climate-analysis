# env variables used in the script for writing csv file
export NUM_NODES=4
export PLACE=scatter:excl
export THREADS
export PROCESSES
export CPUS_PER_NODE
export PROBLEM_SIZE

# load modules
module load mpich-3.2 hdf5-1.10.5--gcc-9.1.0 netcdf-4.7.0--gcc-9.1.0

# make file
mpicc -std=c99 -Wall -g -fopenmp -o ./parallel.out ./parallel.c -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -lm

for PROBLEM_SIZE in 0.25 0.5 1
do
    for PROCESSES in 4 8 16 24 #32 64 96 128
    do
        for THREADS in 1 4
        do
            for RUN in $(seq 1 3)
            do
                mkdir -p output/log/mpi_openMP/problem_${PROBLEM_SIZE}/np_${PROCESSES}/thr_${THREADS}/run_${RUN}
                export CPUS_PER_NODE=$((($THREADS*$PROCESSES)/$NUM_NODES))
                qsub -Wblock=true -N mpi_openmp_run${RUN} -l select=${NUM_NODES}:ncpus=$(($CPUS_PER_NODE)):mem=4gb -l place=${PLACE} -o output/log/mpi_openMP/problem_${PROBLEM_SIZE}/np_${PROCESSES}/thr_${THREADS}/run_${RUN}/result.log -e output/log/mpi_openMP/problem_${PROBLEM_SIZE}/np_${PROCESSES}/thr_${THREADS}/run_${RUN}/error.log parallel.sh
            done
        done
    done
done