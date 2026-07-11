// include headers
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

// start of program
int main(int argc, char* argv[]) {

    // declare variables 
    int rank;                       // value corresponding to this MPI process
    int size;                       // total number of MPI processes
    const int data_size = 16384;    // set data size as a constant (must be larger than internal MPI buffer)
    double data_send[data_size];    // the data the rank sends
    double data_recv[data_size];    // the data the rank recieves
    int partner;                    // the rank of the other process

    // initialise the MPI implementation, determine size and rank 
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // check there are two MPI processes
    if (size != 2) {
        fprintf(stderr, "This program must be launched with 2 MPI processes, instead of %i\n", size);
        return EXIT_FAILURE;
    }

    // set initial value to rank
    for (int i=0; i<data_size; i++) { data_send[i] = (double)rank; }

    // set partner to other rank
    partner = (rank + 1) % 2;

    /* ===================================================================== */
    /* Step 1: Remove the deadlock
     *    Summary:
     *      Reorder the send and recv calls to remove the deadlock.
     *
     *    Hint:
     *      One rank must call MPI_Send when the other is calling MPI_Recv.
     *      
     */
    /* ===================================================================== */ 

    // TODO: which rank sends first?
    if (rank == 0) {

        // send data and print that it has been sent
        MPI_Send(data_send, data_size, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD); 
        printf("Rank %i sent %lf, ...\n", rank, data_send[0]);

        // recieve data and print that it has been received
        MPI_Recv(data_recv, data_size, MPI_DOUBLE, partner, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Rank %i received %lf, ...\n", rank, data_recv[0]);

    // TODO: what does the other rank do to match?
    } else {

        // recieve data and print that it has been received
        MPI_Recv(data_recv, data_size, MPI_DOUBLE, partner, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Rank %i received %lf, ...\n", rank, data_recv[0]);

        // send data and print that it has been sent
        MPI_Send(data_send, data_size, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD); 
        printf("Rank %i sent %lf, ...\n", rank, data_send[0]);

    }

    // finalise the MPI implementation 
    MPI_Finalize();
 
    // exit program
    return EXIT_SUCCESS;

}
