

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

 
fftw_plan r2c_1d_plan, c2r_1d_plan;

	
int main(int argc, char **argv)
{
	const int N0 =8;
	 
	Array<complex<double>,1>  *A;




	A = new Array<complex<double>,1>(N0);
		

	r2c_1d_plan = fftw_plan_dft_r2c_1d(N0, (DP*)((*A).data()), 
											(fftw_complex*)((*A).data()), FFTW_MEASURE);
	
	c2r_1d_plan = fftw_plan_dft_c2r_1d(N0,(fftw_complex*)((*A).data()), 
											(DP*)((*A).data()), FFTW_MEASURE);

	*A = 0.0;

	real((*A)(1)) = 1.5;
	
	cout << *A << endl;
	
	fftw_execute_dft_c2r(c2r_1d_plan, reinterpret_cast<fftw_complex*>((*A).data()), 
		   reinterpret_cast<double*>((*A).data()));
		
	cout << *A << endl;
	 
	fftw_execute_dft_r2c(r2c_1d_plan, reinterpret_cast<double*>((*A).data()), 
			 reinterpret_cast<fftw_complex*>((*A).data()));
 
	
	 cout << "AFTER FFT, A "   << (*A) << endl;
}

