

#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <string>

#ifdef SEEK_SET
#undef SEEK_SET
#endif
#ifdef SEEK_CUR
#undef SEEK_CUR
#endif
#ifdef SEEK_END
#undef SEEK_END
#endif

#include <mpi.h>
#include <fftw3-mpi.h>

#ifndef DP
#define DP double
#endif

#ifndef complx
#define complx  complex<DP>
#endif

// #include <iostream>
// #include <blitz/array.h>
// #include <complex>
// #include <cmath>
// #include <fftw3_mpi.h>
//#include "fourier.h"

using namespace blitz;

int my_id;								// My process id
int numprocs;							// No of processors
const int master_id=0;						// Id of master proc
 int local_N1, local_N1_start;
  
 
fftw_plan r2c_plan_FOUR, c2r_plan_FOUR;
	
int main(int argc, char **argv)
{

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	
	MPI_Status status;
	
	 fftw_mpi_init();
    int N[4] = {0,4,4,4}; 
	  
	/* get local data size and allocate */
	int alloc_local = fftw_mpi_local_size_3d(N[1], N[2], N[3]/2+1, MPI_COMM_WORLD,
                                              &local_N1, &local_N1_start);
											  
	
	Array<complex<double>,3> A(local_N1, N[2], N[3]/2+1);
	Array<complex<double>,3> B(N[2], local_N1,  N[3]/2+1);
	Array<complx,2> temp_plane(N[2], N[3]/2+1);
	
	int data_size = 2* N[2] *(N[3]/2+1);
	
	for (int i=0; i<local_N1; i++)
		real(A(i,0,0)) = (my_id+1)*(i+1);
	
		cout << A << endl;
			
	if (my_id > 0) {												
		temp_plane = A(0, Range::all(),Range::all());
		MPI_Send( reinterpret_cast<double*>(temp_plane.data()), data_size, MPI_DOUBLE, 
						my_id-1, my_id, MPI_COMM_WORLD );	
	}
						
	for (int i = 1; i <= local_N1-1; i++)   { 
		temp_plane = A(i,Range::all(),Range::all());
		A(i-1,Range::all(),Range::all()) = temp_plane;
	}
	
	if (my_id < numprocs-1) {
		MPI_Recv( reinterpret_cast<double*>(temp_plane.data()), data_size, 
						MPI_DOUBLE, my_id+1, my_id+1, MPI_COMM_WORLD, &status );							
		A(local_N1-1,Range::all(),Range::all()) = temp_plane;
	}
	else		// my_id = numproc
		A(local_N1-1,Range::all(),Range::all()) = 0.0;
	
	cout << A << endl;
	
	 B = A.transpose(1,0,2);

	cout << "B" << B << endl;

		 
/*
	string hi, hi2;
	
	char msg[20]; 
	
	 if (my_id == 0) {
		hi = "HELLO S";
	
		strcpy(msg, hi.c_str());
	}

	MPI_Bcast(msg, 20, MPI_CHAR, 0, MPI_COMM_WORLD );
	
	hi2 = msg;
	
	cout << hi2 << " "  << hi2.length() << endl;
	*/

  MPI_Finalize();
}


