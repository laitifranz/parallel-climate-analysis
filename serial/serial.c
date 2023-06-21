// Code adapted from 
//   https://docs.unidata.ucar.edu/netcdf-c/current/pres__temp__4D__rd_8c_source.html for reading task
//   https://docs.unidata.ucar.edu/netcdf-c/current/pres__temp__4D__wr_8c.html for writing task

#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <unistd.h>

//#define LOG /*Print log on command line*/

#define WRITE_CSV

#define INPUT_FILE "../merged_file.nc" /*Change input netCDF path here*/
#define OUTPUT_FILE "output/pr_reduced.nc" /*Define output name file. Saved in the same directory of serial.c*/
#define CSV_FILE "output/performance_benchmarks.csv"

#define NLAT 192
#define NLON 288
#define START_LAT -90.0
#define START_LON 0.0
#define LAT_NAME "lat"
#define LON_NAME "lon"
#define TIME_NAME "time"
#define PR_NAME "pr"
#define DEGREES_EAST "degrees_east"
#define DEGREES_NORTH "degrees_north"
#define LAT_UNITS "degrees_north"
#define LON_UNITS "degrees_east"
#define UNITS "units"
#define PR_UNITS "kg m-2 s-1"

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

int main(){
    
    MPI_Init(NULL, NULL);
    FILE *fpt;

    #ifdef WRITE_CSV
    if(access(CSV_FILE, F_OK) == -1){
        fpt = fopen(CSV_FILE, "a");
        fprintf(fpt, "problem_size,time\n");
    }
    else fpt = fopen(CSV_FILE, "a");
    #endif

    double start_time, end_time;
    int retval;

    /* netCDF file ID and variable ID */
    int ncid;
    int pr_varid, lat_varid, lon_varid;
    int lon_dimid, lat_dimid, time_dimid;
    int ndims;
    size_t nrecord_fromID;
    
    /* --- GET INFO STEP --- */
    /* Open the file with read-only access, indicated by NC_NOWRITE flag */
    start_time = MPI_Wtime();

    if ((retval = nc_open(INPUT_FILE, NC_NOWRITE, &ncid))) ERR(retval);
    /* Retrieve dimensions info */
    if ((retval = nc_inq_ndims(ncid, &ndims))) ERR(retval);

    if ((retval = nc_inq_dimid(ncid, LAT_NAME, &lat_dimid))) ERR(retval);
    if ((retval = nc_inq_dimid(ncid, LON_NAME, &lon_dimid))) ERR(retval);
    if ((retval = nc_inq_dimid(ncid, TIME_NAME, &time_dimid))) ERR(retval);

    if ((retval = nc_inq_dimlen(ncid, time_dimid, &nrecord_fromID))) ERR(retval);

    printf("--- INFO: found dim = %d and nrecord = %zu ---\n", ndims, nrecord_fromID);
    
    
    int nrecord = (int) nrecord_fromID;
    float problem_size;
    if(getenv("PROBLEM_SIZE")){
        problem_size = atof(getenv("PROBLEM_SIZE"));
        nrecord = nrecord * problem_size;
        printf("--- INFO: new problem size set: nrecord = %d with problem size = %f ---\n", nrecord, problem_size);
    }
    ndims -= 1; // we do not use the bnds dimension, so we need to reduce by one the total number of dims

    /* Retrieve variables info */
    if ((retval = nc_inq_varid(ncid, LAT_NAME, &lat_varid))) ERR(retval);
    if ((retval = nc_inq_varid(ncid, LON_NAME, &lon_varid))) ERR(retval);
    if ((retval = nc_inq_varid(ncid, PR_NAME, &pr_varid))) ERR(retval);

    float lats[NLAT], lons[NLON];
    if ((retval = nc_get_var_float(ncid, lat_varid, &lats[0]))) ERR(retval);  /* Read the coordinate variable data. */
    if ((retval = nc_get_var_float(ncid, lon_varid, &lons[0]))) ERR(retval);
    
    end_time = MPI_Wtime();
    printf("## GET INFO STEP: %f seconds ##\n", end_time - start_time);
    
    #ifdef LOG
    /* DEBUG: Check if lats and lons were read correctly */
    for (int lat = 0; lat < NLAT; lat++)
        printf("%f \n", lats[lat]);
    for (int lon = 0; lon < NLON; lon++)
        printf("%f \n", lons[lon]);
    #endif


    /* --- READING STEP --- */
    start_time = MPI_Wtime();

    float pr_in[NLAT][NLON], pr_out[NLAT][NLON];
    size_t start[] = {0, 0, 0}; // 3-dim array
    size_t count[] = {1, NLAT, NLON}; // 3-dim array 

    /* Get the precipitation value for each time step for each point of the grid, then sum it up to obtain the sum over the years */
    for (int rec = 0; rec < nrecord; rec++){
        start[0] = rec;
        if ((retval = nc_get_vara_float(ncid, pr_varid, start, count, &pr_in[0][0]))) ERR(retval);

        for(int lat = 0; lat < NLAT; lat++){
            for(int lon = 0; lon < NLON; lon++){
                pr_out[lat][lon] = pr_out[lat][lon] + pr_in[lat][lon];
            }
        }
    }

    end_time = MPI_Wtime();

    #ifdef WRITE_CSV
    fprintf(fpt, "%f, %f\n", problem_size, end_time - start_time);
    #endif

    /* Get the average precipitations over the years for each point of the grid (average precipitations) */
    for(int lat = 0; lat < NLAT; lat++){
        for(int lon = 0; lon < NLON ; lon++){
            pr_out[lat][lon] /= nrecord;
        }
    }

    if ((retval = nc_close(ncid))) ERR(retval);     /* Close the file */
    printf("--- SUCCESS reading data from file %s ---\n", INPUT_FILE);


    /* --- WRITING STEP --- */

    start_time = MPI_Wtime();

    /* Create the file */
    if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid))) ERR(retval);

    if ((retval = nc_def_dim(ncid, LAT_NAME, NLAT, &lat_dimid))) ERR(retval);
    if ((retval = nc_def_dim(ncid, LON_NAME, NLON, &lon_dimid))) ERR(retval);

    if ((retval = nc_def_var(ncid, LAT_NAME, NC_FLOAT, 1, &lat_dimid, &lat_varid))) ERR(retval);
    if ((retval = nc_def_var(ncid, LON_NAME, NC_FLOAT, 1, &lon_dimid, &lon_varid))) ERR(retval);
    int dimids[] = {lat_dimid, lon_dimid}; // 2-dim array
    ndims -= 1; // we do not use the time dimension, so we need to reduce by one the total number of dims
    if ((retval = nc_def_var(ncid, PR_NAME, NC_FLOAT, ndims, dimids, &pr_varid))) ERR(retval); /* Define the netCDF variables for the precipitation data */

    if ((retval = nc_put_att_text(ncid, lat_varid, UNITS, strlen(DEGREES_NORTH), DEGREES_NORTH))) ERR(retval);
    if ((retval = nc_put_att_text(ncid, lon_varid, UNITS, strlen(DEGREES_EAST), DEGREES_EAST))) ERR(retval);
    if ((retval = nc_put_att_text(ncid, pr_varid, UNITS, strlen(PR_UNITS), PR_UNITS))) ERR(retval);     /* Assign units attributes to the netCDF variables */

    if ((retval = nc_enddef(ncid))) ERR(retval); /* End define mode */

    /* Write the coordinate variable data. This will put the latitudes and longitudes of our data grid into the netCDF file */
    if ((retval = nc_put_var_float(ncid, lat_varid, &lats[0]))) ERR(retval);
    if ((retval = nc_put_var_float(ncid, lon_varid, &lons[0]))) ERR(retval);

    size_t count_write[] = {NLAT, NLON}; // 2-dim array
    size_t start_write[] = {0, 0}; // 2-dim array 

    /* Write in the new netCDF file the average over the years for each point of the grid of the precipitations */
    if ((retval = nc_put_vara_float(ncid, pr_varid, start_write, count_write, &pr_out[0][0]))) ERR(retval);

    end_time = MPI_Wtime();

    /* DEBUG: time for execution*/
    printf("## WRITING STEP: %f seconds ##\n", end_time - start_time);

    /*Close the file, freeing all resources */
    if ((retval = nc_close(ncid))) ERR(retval);       

    printf("--- SUCCESS writing data on file %s ---\n", OUTPUT_FILE);

    MPI_Finalize();
    return 0;
}