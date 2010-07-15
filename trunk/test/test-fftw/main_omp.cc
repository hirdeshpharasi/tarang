

#include <iostream>
// #include <blitz/array.h>
// #include <complex>
// #include <cmath>

#include "fourier.h"

fftw_plan r2c_plan_FOUR, c2r_plan_FOUR;

	
main()
{
  int dim=2;
  int N[3]={0,8,8};
  
  Array<complx,2> A(N[1], N[2]/2+1);           // Temp location to hold *Vi(:,j)

  A(0,1)=complex<DP>(1.0,0); 
    
  A(0,2)=complex<DP>(0.3,0); // (1,0,1) 1/4i
 
  
  
	fftw_init_threads();
	fftw_plan_with_nthreads(2); 
	
	Init_fftw_plan_FOUR(N, A);
	
	for (int i=1; i<=10; i++) {
	ArrayFFT_FOUR(c2r_plan_FOUR, N, A);
	cout << A << endl;
	
	ArrayIFFT_FOUR(r2c_plan_FOUR, N, A);
	cout << A << endl;
	}
	
	fftw_cleanup_threads();

} 



