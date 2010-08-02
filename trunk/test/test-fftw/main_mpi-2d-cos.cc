

#include <blitz/array.h>
#include <complex>
#include <cmath>

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
ptrdiff_t local_n0, local_0_start;
  
 
fftw_plan r2c_plan_FOUR, c2r_plan_FOUR, plan;

DP Local_abs_sqr_FOUR(int N[], Array<complx,3> A);
DP Total_abs_sqr_FOUR(int N[], Array<complx,3> A);
	
int main(int argc, char **argv)
{
	const int N0 =8, N1 = 8;
	int NN[3]; 

	ptrdiff_t alloc_local;
	int i, j;
	DP tot;
	 
	Array<complex<double>,2>  *A;
	Array<complex<double>,2>  *B;
 
	NN[1] = N0; NN[2] = N1; 

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	 
	cout << "MPI " << my_id << " "  << numprocs << endl;

	fftw_mpi_init();
 
	 /* get local data size and allocate */
	alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
										  &local_n0, &local_0_start);
 
	cout << "Myid , Localn0, localnostart:   " << my_id << " " << local_n0 <<" "  << local_0_start << endl;
	
	/*

	A = new Array<complex<double>,2>(local_n0, N1);

	B = new Array<complex<double>,2>(local_n0, N1);
	

	*A = 0.0;
	if (my_id == 0)
	{
		real((*A)(2,2)) = 5.0;  imag((*A)(2,2)) = 5.0;
	}
	
		 
	cout << my_id << " " << *A << endl;
		
	plan = fftw_mpi_plan_dft_2d(N0, N1,  reinterpret_cast<fftw_complex*>((*A).data()), 
									reinterpret_cast<fftw_complex*>((*A).data()), MPI_COMM_WORLD,
									FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftw_execute(plan);
		 
		 cout << my_id << " " << *A << endl;
	
	*/
	

	A = new Array<complex<double>,2>(local_n0, N1/2+1);
	
	B = new Array<complex<double>,2>(local_n0, N1/2+1);
	

	r2c_plan_FOUR = fftw_mpi_plan_dft_r2c_2d(NN[1], NN[2], 
					reinterpret_cast<double*>((*A).data()), 
					reinterpret_cast<fftw_complex*>((*A).data()),  
					MPI_COMM_WORLD, FFTW_MEASURE);
	
	c2r_plan_FOUR = fftw_mpi_plan_dft_c2r_2d(NN[1], NN[2],
				reinterpret_cast<fftw_complex*>((*A).data()), 
				reinterpret_cast<double*>((*A).data()), MPI_COMM_WORLD, FFTW_MEASURE);

	*A = 0.0;
	if (my_id == 0)
	{
		real((*A)(1,1)) = 1.0;  imag((*A)(1,1)) = 1.0;
	} 
	
	cout << my_id << " " << *A << endl;
	
	fftw_execute_dft_c2r(c2r_plan_FOUR, reinterpret_cast<fftw_complex*>((*A).data()), 
		   reinterpret_cast<double*>((*A).data()));
		
	cout << my_id << " " << *A << endl;
	 
	fftw_execute_dft_r2c(r2c_plan_FOUR, reinterpret_cast<double*>((*A).data()), 
			 reinterpret_cast<fftw_complex*>((*A).data()));
	
	cout << my_id << " " << *A << endl;
 
	
	 cout << "AFTER FFT, myid, A " << my_id <<  *A << endl;
	 
	 
	 fftw_destroy_plan(r2c_plan_FOUR);
	 fftw_destroy_plan(c2r_plan_FOUR);	
	 MPI_Finalize();
}

