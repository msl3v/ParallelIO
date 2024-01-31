/*
 * Tests for PIO distributed arrays. This test uses 1 dimension,
 * everything very simple. ;-)
 *
 * @author Ed Hartnett
 * @date 2/27/17
 */
#include <config.h>
#include <pio.h>
#include <pio_internal.h>
#include <pio_tests.h>

/* The number of tasks this test should run on. */
#define TARGET_NTASKS 4

/* The minimum number of tasks this test should run on. */
#define MIN_NTASKS 4

/* The name of this test. */
#define TEST_NAME "test_darray_1d"

/* Number of processors that will do IO. */
#define NUM_IO_PROCS 1

/* Number of computational components to create. */
#define COMPONENT_COUNT 1

/* The number of dimensions in the example data. In this test, we
 * are using three-dimensional data. */
#define NDIM 1

/* The length of our sample data along the dimension. */
#define DIM_LEN 4

/* This is the length of the map for each task. */
#define EXPECTED_MAPLEN 2

/* The number of timesteps of data to write. */
#define NUM_TIMESTEPS 2

/* The name of the variable in the netCDF output files. */
#define VAR_NAME "RedShirtSurvival"

/* The dimension names. */
#define DIM_NAME "episode"
#define DIM_NAME_2 "phaser_draws"

/* Create a 1D decomposition.
 *
 * @param ntasks the number of available tasks
 * @param my_rank rank of this task.
 * @param iosysid the IO system ID.
 * @param dim_len an array of length 3 with the dimension sizes.
 * @param ioid a pointer that gets the ID of this decomposition.
 * @param pio_type the type that will be used for basetype.
 * @returns 0 for success, error code otherwise.
 **/
int create_decomposition_1d(int ntasks, int my_rank, int iosysid, int pio_type, int *ioid)
{
    PIO_Offset elements_per_pe;     /* Array elements per processing unit. */
    int dim_len_1d[NDIM] = {DIM_LEN};
    int ret;

    /* How many data elements per task? In this example we will end up
     * with 2. */
    elements_per_pe = DIM_LEN / ntasks;

    PIO_Offset compdof[elements_per_pe];

    /* Don't forget to add 1! */
    compdof[0] = my_rank + 1;

    /* This means fill value will be used here. */
    compdof[1] = 0;

    /* Create the PIO decomposition for this test. */
    if ((ret = PIOc_InitDecomp(iosysid, pio_type, NDIM, dim_len_1d, elements_per_pe,
                               compdof, ioid, NULL, NULL, NULL)))
        ERR(ret);

    return 0;
}

/**
 * Test fill values and darrays.
 *
 * @param iosysid the IO system ID.
 * @param ioid the ID of the decomposition.
 * @param pio_type the type of the data.
 * @param num_flavors the number of IOTYPES available in this build.
 * @param flavor array of available iotypes.
 * @param my_rank rank of this task.
 * @param test_comm the MPI communicator running the test.
 * @returns 0 for success, error code otherwise.
 */
int test_darray_fill(int iosysid, int ioid, int pio_type, int num_flavors, int *flavor,
                     int my_rank, MPI_Comm test_comm)
{
#define NUM_FILLVALUE_PRESENT_TESTS 2
    char filename[PIO_MAX_NAME + 1]; /* Name for the output files. */
    int dimid;     /* The dimension ID. */
    int ncid;      /* The ncid of the netCDF file. */
    int varid;     /* The ID of the netCDF varable. */
    PIO_Offset arraylen = 2;
    void *test_data;
    void *fillvalue;
    void *test_data_in;
    void *expected_in;
    PIO_Offset type_size;             /* Size of the data type. */
    /* My rank as each type. */
    signed char my_byte_rank = my_rank;
    char my_char_rank = my_rank;
    short my_short_rank = my_rank;
    float my_float_rank = my_rank;
    double my_double_rank = my_rank;
#ifdef _NETCDF4
    unsigned char my_ubyte_rank = my_rank;
    unsigned short my_ushort_rank = my_rank;
    unsigned int my_uint_rank = my_rank;
    long long my_int64_rank = my_rank;
    unsigned long long my_uint64_rank = my_rank;
#endif /* _NETCDF4 */

    /* Default fill value for each type. */
    signed char byte_fill = NC_FILL_BYTE;
    char char_fill = NC_FILL_CHAR;
    short short_fill = NC_FILL_SHORT;
    int int_fill = NC_FILL_INT;
    float float_fill = NC_FILL_FLOAT;
    double double_fill = NC_FILL_DOUBLE;
#ifdef _NETCDF4
    unsigned char ubyte_fill = NC_FILL_UBYTE;
    unsigned short ushort_fill = NC_FILL_USHORT;
    unsigned int uint_fill = NC_FILL_UINT;
    long long int64_fill = NC_FILL_INT64;
    unsigned long long uint64_fill = NC_FILL_UINT64;
#endif /* _NETCDF4 */

    void *bufr;
    int ret; /* Return code. */

    GDALDatasetH hDSp;
    int iotype = PIO_IOTYPE_GDAL;

    /* Use PIO to create the example file in each of the four
     * available ways. */
    for (int fmt = 2; fmt < 3; fmt++)
    {
      flavor[fmt] = 2;
        /* BYTE and CHAR don't work with pnetcdf. Don't know why yet. */
/*        if (flavor[fmt] == PIO_IOTYPE_PNETCDF && (pio_type == PIO_BYTE || pio_type == PIO_CHAR))
            continue;
*/
        /* NetCDF-4 types only work with netCDF-4 formats. */
        if (pio_type > PIO_DOUBLE && flavor[fmt] != PIO_IOTYPE_NETCDF4C &&
            flavor[fmt] != PIO_IOTYPE_NETCDF4P)
            continue;

	/* Create the filename. */
	sprintf(filename, "data/cb_2018_us_region_20m.shp");

	/* Reopen the file. */
	printf("Here1\n");
	if ((ret = GDALc_openfile(iosysid, &ncid, &hDSp, &iotype, filename, PIO_NOWRITE))) {
	  printf("Error %d\n",ret);
	  ERR(ret);
	}

	printf("Here2\n");
	if ((ret = GDALc_inq_fieldid(ncid, "GEOID", &varid)))
	  ERR(ret);

	printf("Here3\n");
	/* Allocate space for data. */
	if (!(test_data_in = malloc(sizeof(double) * arraylen)))
	  ERR(PIO_ENOMEM);

	/* Read the data. */
	if ((ret = PIOc_read_darray(ncid, varid, ioid, arraylen, test_data_in)))
	  ERR(ret);

	printf("test_data_in[%d] %f\n",my_rank,((double *)test_data_in)[0]);

	/* Check the (first) result. */
	//            if (memcmp(test_data_in, expected_in, type_size))
	//                return ERR_WRONG;

	/* Free resources. */
	free(test_data_in);

	/* Close the netCDF file. */
	if ((ret = PIOc_closefile(ncid)))
	  ERR(ret);
    } /* next iotype */

    return PIO_NOERR;
}

/* Run tests for darray functions. */
int main(int argc, char **argv)
{
#define NUM_REARRANGERS_TO_TEST 1
    int rearranger[NUM_REARRANGERS_TO_TEST] = {PIO_REARR_BOX};//, PIO_REARR_SUBSET};
//#ifdef _NETCDF4
//#define NUM_TYPES_TO_TEST 11
//    int test_type[NUM_TYPES_TO_TEST] = {PIO_BYTE, PIO_CHAR, PIO_SHORT, PIO_INT, PIO_FLOAT, PIO_DOUBLE,
//                                        PIO_UBYTE, PIO_USHORT, PIO_UINT, PIO_INT64, PIO_UINT64};
//#else
#define NUM_TYPES_TO_TEST 1
    int test_type[NUM_TYPES_TO_TEST] = {PIO_DOUBLE};//{PIO_BYTE, PIO_CHAR, PIO_SHORT, PIO_INT, PIO_FLOAT, PIO_DOUBLE};
//#endif /* _NETCDF4 */
    int my_rank;
    int ntasks;
    int num_flavors; /* Number of PIO netCDF flavors in this build. */
    int flavor[NUM_FLAVORS]; /* iotypes for the supported netCDF IO flavors. */
    MPI_Comm test_comm; /* A communicator for this test. */
    int ret;         /* Return code. */

    /* Initialize test. */
    if ((ret = pio_test_init2(argc, argv, &my_rank, &ntasks, MIN_NTASKS,
                              MIN_NTASKS, -1, &test_comm)))
        ERR(ERR_INIT);

    OGRRegisterAll();
    PIOc_set_log_level(4);

    if ((ret = PIOc_set_iosystem_error_handling(PIO_DEFAULT, PIO_RETURN_ERROR, NULL)))
        return ret;

    /* Only do something on max_ntasks tasks. */
    if (my_rank < TARGET_NTASKS)
    {
        int iosysid;  /* The ID for the parallel I/O system. */
        int ioid;     /* Decomposition ID. */
        int ioproc_stride = 1;    /* Stride in the mpi rank between io tasks. */
        int ioproc_start = 0;     /* Rank of first processor to be used for I/O. */
        int ret;      /* Return code. */

        /* Figure out iotypes. */
        if ((ret = get_iotypes(&num_flavors, flavor)))
            ERR(ret);

        for (int r = 0; r < NUM_REARRANGERS_TO_TEST; r++)
        {
            /* Initialize the PIO IO system. This specifies how many and
             * which processors are involved in I/O. */
            if ((ret = PIOc_Init_Intracomm(test_comm, TARGET_NTASKS, ioproc_stride,
                                           ioproc_start, rearranger[r], &iosysid)))
                return ret;

            /* Run tests for each data type. */
            for (int t = 0; t < NUM_TYPES_TO_TEST; t++)
            {
                /* Decompose the data over the tasks. */
                if ((ret = create_decomposition_1d(TARGET_NTASKS, my_rank, iosysid, test_type[t],
                                                   &ioid)))
                    return ret;

                /* Run tests. */
                if ((ret = test_darray_fill(iosysid, ioid, test_type[t], num_flavors, flavor,
                                            my_rank, test_comm)))
                    return ret;

                /* Free the PIO decomposition. */
                if ((ret = PIOc_freedecomp(iosysid, ioid)))
                    ERR(ret);
            }

            /* Finalize PIO system. */
            if ((ret = PIOc_free_iosystem(iosysid)))
                return ret;
        } /* next rearranger */

    } /* endif my_rank < TARGET_NTASKS */

    /* Finalize the MPI library. */
    if ((ret = pio_test_finalize(&test_comm)))
        return ret;

    printf("%d %s SUCCESS!!\n", my_rank, TEST_NAME);
    return 0;
}
