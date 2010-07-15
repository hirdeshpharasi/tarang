
/* Uses fftw3; Real_to_complex transform.; blitz++
   Before FFT, zero_pad.
   After FFT, normalize (divide by N1*N2*...)
*/

#include "fft_my.h"

using namespace blitz ;

extern fftw_plan r2c, c2r;

// init_fftw
void init_fftw(int NN[], Array<complx,1> A)
{
  r2c = fftw_plan_dft_r2c_1d(NN[1],(double*)(A.data()), 
			     (fftw_complex*)(A.data()),  FFTW_ESTIMATE);
  c2r = fftw_plan_dft_c2r_1d(NN[1],(fftw_complex*)(A.data()), 
			     (double*)(A.data()),FFTW_ESTIMATE);
}

void init_fftw(int NN[], Array<complx,2> A)
{
  r2c = fftw_plan_dft_r2c_2d(NN[1],NN[2],(double*)(A.data()), 
			     (fftw_complex*)(A.data()),  FFTW_ESTIMATE);
  c2r = fftw_plan_dft_c2r_2d(NN[1],NN[2],(fftw_complex*)(A.data()), 
			     (double*)(A.data()),FFTW_ESTIMATE);
}

void init_fftw(int NN[], Array<complx,3> A)
{
  r2c = fftw_plan_dft_r2c_3d(NN[1],NN[2],NN[3],(double*)(A.data()), 
			     (fftw_complex*)(A.data()),  FFTW_ESTIMATE);
  c2r = fftw_plan_dft_c2r_3d(NN[1],NN[2],NN[3],(fftw_complex*)(A.data()), 
			     (double*)(A.data()),FFTW_ESTIMATE);
			     
			    cout << "in fft_int " << endl;
}

// FFT
void arrayFFT(fftw_plan r2c, Array<complx,1> A)
{
  fftw_execute_dft_r2c(r2c, reinterpret_cast<double*>(A.data()), 
			     reinterpret_cast<fftw_complex*>(A.data()));
}
void arrayFFT(fftw_plan r2c, Array<complx,2> A)
{
  fftw_execute_dft_r2c(r2c, reinterpret_cast<double*>(A.data()), 
		       reinterpret_cast<fftw_complex*>(A.data()));
}
void arrayFFT(fftw_plan r2c, Array<complx,3> A) 
{
  fftw_execute_dft_r2c(r2c, reinterpret_cast<double*>(A.data()), 
		       reinterpret_cast<fftw_complex*>(A.data()));
}

//IFFT
void arrayIFFT(fftw_plan c2r, Array<complx,1> A)
{
  fftw_execute_dft_c2r(c2r, reinterpret_cast<fftw_complex*>(A.data()), 
		       reinterpret_cast<double*>(A.data()));
}

void arrayIFFT(fftw_plan c2r, Array<complx,2> A)
{
  fftw_execute_dft_c2r(c2r, reinterpret_cast<fftw_complex*>(A.data()), 
		       reinterpret_cast<double*>(A.data()));
}

void arrayIFFT(fftw_plan c2r, Array<complx,3> A)
{
  fftw_execute_dft_c2r(c2r, reinterpret_cast<fftw_complex*>(A.data()), 
		       reinterpret_cast<double*>(A.data()));
}


