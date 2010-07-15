
/**************************** fft_my.h *********************/
// fft_my.h

#include <blitz/array.h>
#include <complex>
#include <fftw3.h>

#ifndef DP
#define DP double
#endif

#ifndef complx
#define complx  complex<DP>
#endif

using namespace blitz ;

void init_fftw(int NN[], Array<complx,1> A);
void init_fftw(int NN[], Array<complx,2> A);
void init_fftw(int NN[], Array<complx,3> A);

void arrayFFT(fftw_plan ,Array<complx,1> A);
void arrayFFT(fftw_plan ,Array<complx,2> A);
void arrayFFT(fftw_plan ,Array<complx,3> A); 

void arrayIFFT(fftw_plan, Array<complx,1> A);
void arrayIFFT(fftw_plan, Array<complx,2> A);
void arrayIFFT(fftw_plan, Array<complx,3> A);

void fft_norm(int N[], Array<complx,1> A);
void fft_norm(int N[], Array<complx,2> A);
void fft_norm(int N[], Array<complx,3> A);


/**************************** fft_my.h *********************/
