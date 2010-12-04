#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <complex.h>
#include <mpi.h>
#include <fftw3-mpi.h>
#include <time.h>

using namespace std;

int main(int argc, char** argv)
{
    const int size_of_array = 64; const int dim = 3; const int number_of_transforms = 10;
    unsigned int numBytes;
    int my_id;
    int numprocs;
    const int master_id = 0;
    ptrdiff_t local_n0, local_0_start, alloc_local;
    time_t start, end;

    fftw_complex *data_array;
    fftw_plan forward_plan, backward_plan;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    fftw_mpi_init();



    switch (dim)
    {
        case 2:
            alloc_local = fftw_mpi_local_size_2d(size_of_array, size_of_array, MPI_COMM_WORLD, &local_n0, &local_0_start);
            numBytes = pow((double)size_of_array, (double)dim)*sizeof(fftw_complex)*alloc_local;
            data_array = (fftw_complex*) fftw_malloc(numBytes);
            forward_plan = fftw_mpi_plan_dft_r2c_2d(size_of_array, size_of_array, (double*)data_array, data_array, MPI_COMM_WORLD, FFTW_PATIENT);
            backward_plan = fftw_mpi_plan_dft_c2r_2d(size_of_array, size_of_array, data_array, (double*)data_array, MPI_COMM_WORLD, FFTW_PATIENT);

            for (int i=0; i < local_n0; i++)
            {
                for (int j=0; j < size_of_array; j++)
                {
//                    data_array[j + size_of_array*i] = rand()/(double)RAND_MAX;
                }
            }
            break;

        case 3:
            alloc_local = fftw_mpi_local_size_3d(size_of_array, size_of_array, size_of_array, MPI_COMM_WORLD, &local_n0, &local_0_start);
            numBytes = pow((double)size_of_array, (double)dim)*sizeof(fftw_complex)*alloc_local;
            data_array = (fftw_complex*) fftw_malloc(numBytes);
            forward_plan = fftw_mpi_plan_dft_r2c_3d(size_of_array, size_of_array, size_of_array, (double*)data_array, data_array, MPI_COMM_WORLD, FFTW_PATIENT);
            backward_plan = fftw_mpi_plan_dft_c2r_3d(size_of_array, size_of_array, size_of_array, data_array, (double*)data_array, MPI_COMM_WORLD, FFTW_PATIENT);

            for (int i=0; i < local_n0; i++)
            {
                for (int j=0; j < size_of_array; j++)
                {
                    for (int k=0; k < size_of_array; k++)
                    {
                        //data_array[k + size_of_array*(j + size_of_array*i)] = rand()/(double)RAND_MAX;
                    }
                }
            }
            break;
    }

    time(&start);

    for (int i = 0; i < number_of_transforms; i++) 
    {
        fftw_execute(forward_plan); 
        fftw_execute(backward_plan); 
    }

    time(&end);

    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);
    fftw_free(data_array);

    MPI_Finalize();

//    if (my_id==master_id)
//    {
//        cout << endl << "PROGRAM TERMINATING HERE" << endl << endl·
//        << "TOTAL TIME ELAPSED: " << difftime(end, start) << " sec" << endl;
//    }
}
