export PROBLEM_SIZE

# load modules
module load mpich-3.2 hdf5-1.10.5--gcc-9.1.0 netcdf-4.7.0--gcc-9.1.0

# make file
mpicc -std=c99 -Wall -g -o ./serial.out ./serial.c -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -lm

for PROBLEM_SIZE in 0.25 0.5 1
do
    for RUN in $(seq 1 3)
    do 
        # create directory if not exists
        mkdir -p output/log/problem_${PROBLEM_SIZE}/run_${RUN}

        #launch qsub command
        qsub -Wblock=true -N serial_run${RUN} -o output/log/problem_${PROBLEM_SIZE}/run_${RUN}/result.log -e output/log/problem_${PROBLEM_SIZE}/run_${RUN}/error.log serial.sh
    done
done