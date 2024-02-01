#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>

float **alloc_2d_float( int ndim1, int ndim2 ) {
    float **array2 = malloc( ndim1 * sizeof( float * ) );
    int i;

    if( array2 != NULL ){
        array2[0] = malloc( ndim1 * ndim2 * sizeof( float ) );
        if( array2[ 0 ] != NULL ) {
            for( i = 1; i < ndim1; i++ )
                array2[i] = array2[0] + i * ndim2;
        }

        else {
            free( array2 );
            array2 = NULL;
        }
    }
    return array2;
}

void free_2d_float( float **array ) {
    if (array != NULL) {
        free(array[0]);
        free(array);
    }
    return;
}

void init_array2d(float **array, int ndim1, int ndim2, float data) {
    for (int i=0; i<ndim1; i++) 
        for (int j=0; j<ndim2; j++)
            array[i][j] = data;
    return;
}

void print_array2d(float **array, int ndim1, int ndim2) {
    for (int i=0; i<ndim1; i++) {
        for (int j=0; j<ndim2; j++) {
            printf("%6.2f ", array[i][j]);
        }
        printf("\n");
    }
    return;
}


int main(int argc, char ** argv)
{
    int size, rank;
    int dim_size[2];
    int periods[2];
    MPI_Datatype block_type, resized_type;

    float **array;
    float **whole_array;
    float *recvptr;

    int *counts, *disps;

    int n = 16;
    int rows_per_core;
    int cols_per_core;
    int i, j;

    int whole_array_size[2];
    int sub_array_size[2];
    int starts[2];
    int A, B;

    /* Initialise MPI */
    MPI_Init(&argc, &argv);

    /* Get the rank for this process, and the number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        /* If we're the master process */
        whole_array = alloc_2d_float(n, n);
        recvptr = &(whole_array[0][0]);

        /* Initialise whole array to silly values */
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                whole_array[i][j] = 9999.99;
            }
        }

        print_array2d(whole_array, n, n);
        puts("\n\n");
    }

    /* Create the cartesian communicator */
    MPI_Dims_create(size, 2, dim_size);
    A = dim_size[1];
    B = dim_size[0];
    periods[0] = 1;
    periods[1] = 1;

    rows_per_core = ceil(n / (float) A);
    cols_per_core = ceil(n / (float) B);
    if (rows_per_core*A != n) {
        if (rank == 0) fprintf(stderr,"Aborting: rows %d don't divide by %d evenly\n", n, A);
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    if (cols_per_core*B != n) {
        if (rank == 0) fprintf(stderr,"Aborting: cols %d don't divide by %d evenly\n", n, B);
        MPI_Abort(MPI_COMM_WORLD,2);
    }

    array = alloc_2d_float(rows_per_core, cols_per_core);
    printf("%d, RpC: %d, CpC: %d\n", rank, rows_per_core, cols_per_core);

    whole_array_size[0] = n;             
    sub_array_size  [0] = rows_per_core; 
    whole_array_size[1] = n;
    sub_array_size  [1] = cols_per_core;
    starts[0] = 0; starts[1] = 0;

    MPI_Type_create_subarray(2, whole_array_size, sub_array_size, starts, 
                             MPI_ORDER_C, MPI_FLOAT, &block_type);
    MPI_Type_commit(&block_type);
    MPI_Type_create_resized(block_type, 0, 1*sizeof(float), &resized_type);
    MPI_Type_commit(&resized_type);

    if (array == NULL)
    {
        printf("Problem with array allocation.\nExiting\n");
        MPI_Abort(MPI_COMM_WORLD,3);
    }

    init_array2d(array,rows_per_core,cols_per_core,(float)rank);

    counts = (int *)malloc(size * sizeof(int));
    disps  = (int *)malloc(size * sizeof(int));
    /* note -- we're just using MPI_COMM_WORLD rank here to
     * determine location, not the cart_comm for now... */
    for (int i=0; i<size; i++) {
        counts[i] = 1;   /* one block_type per rank */

        int row = (i % A);
        int col = (i / A);
        /* displacement into the whole_array */
        disps[i] = (col*cols_per_core + row*(rows_per_core)*n);
    }

    MPI_Gatherv(array[0], rows_per_core*cols_per_core, MPI_FLOAT, 
                recvptr, counts, disps, resized_type, 0, MPI_COMM_WORLD);

    free_2d_float(array);
    if (rank == 0) print_array2d(whole_array, n, n);
    if (rank == 0) free_2d_float(whole_array);
    MPI_Finalize();
}