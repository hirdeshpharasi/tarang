

#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <fftw3.h>

#ifdef SEEK_SET
#undef SEEK_SET
#endif
#ifdef SEEK_CUR
#undef SEEK_CUR
#endif
#ifdef SEEK_END
#undef SEEK_END
#endif

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
  


DP Local_abs_sqr_FOUR(int N[], Array<complx,3> A);
DP Total_abs_sqr_FOUR(int N[], Array<complx,3> A);
	
int main(int argc, char **argv)
{
	const int N0 =4, N1 = 3;
	int NN[3]; 

	 
	Array<complex<double>,2>  A(N0, N1);
	Array<complex<double>,2>  B(N0, N1);
 
	NN[1] = N0; NN[2] = N1; 
	
	fftw_plan c2c_2d_forward_plan_FOUR, c2c_2d_inverse_plan_FOUR;
	

	

/*	r2c_plan_FOUR = fftw_mpi_plan_dft_r2c_2d(NN[1], NN[2], 
					reinterpret_cast<double*>((*A).data()), 
					reinterpret_cast<fftw_complex*>((*A).data()),  
					MPI_COMM_WORLD, FFTW_MEASURE);
	
	c2r_plan_FOUR = fftw_mpi_plan_dft_c2r_2d(NN[1], NN[2],
				reinterpret_cast<fftw_complex*>((*A).data()), 
				reinterpret_cast<double*>((*A).data()), MPI_COMM_WORLD, FFTW_MEASURE); */
	
	
	c2c_2d_forward_plan_FOUR = fftw_plan_dft_2d(NN[1], NN[2],
							reinterpret_cast<fftw_complex*>((A).data()), 
							reinterpret_cast<fftw_complex*>((A).data()),	
							FFTW_FORWARD, FFTW_MEASURE); 

	
	c2c_2d_inverse_plan_FOUR = fftw_plan_dft_2d(NN[1], NN[2],
							reinterpret_cast<fftw_complex*>((A).data()), 
							reinterpret_cast<fftw_complex*>((A).data()),		
							FFTW_BACKWARD, FFTW_MEASURE);
	
	A = 0.0;
	
//	real(A(2,2)) = 5.0;  imag(A(2,2)) = 5.0;
	real(A(1,1)) = 5.0; imag(A(1,1)) = 4.0;
	
	cout  << A << endl;
	
//	fftw_execute_dft(c2c_2d_forward_plan_FOUR, reinterpret_cast<fftw_complex*>(A.data()), 
//					 reinterpret_cast<fftw_complex*>(A.data()));
	
	fftw_execute_dft(c2c_2d_inverse_plan_FOUR, reinterpret_cast<fftw_complex*>(A.data()),
					 reinterpret_cast<fftw_complex*>(A.data()));
	
	cout  << A << endl;
	 
}

