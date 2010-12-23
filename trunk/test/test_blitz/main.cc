
#include <iostream> 
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <complex>
#include <cmath>

#ifdef BZ_HAVE_STD
#include <fstream>
#else
#include <fstream.h>
#endif

BZ_USING_NAMESPACE(blitz)

using namespace blitz;

// fftw_plan c2r,r2c;
void xderiv(int NN[],Array<complex<double>,2> A, Array<complex<double>,2> B);

inline double my_pow(double x, int n) {		// n=0,2,4 only
	if (n==0) return 1.0;
	else if (n==2) return (x*x);
	else if (n==4) return (x*x*x*x);
		}


int main()
{
  
  ofstream global_file;	
  global_file.open("test.d");
  
  char* test;
  char* test1;
  
//  test = "HI_S";
// test1 = test;
	cout << test << " " << test1 << endl;
  int dim=3;
  int NN[3]={0,8,8};
  
  Array<complex<double>,2> A(8,5),B(8,5);
  Array<int,3> C(3,4,5);
 //  Array<complex<double>,1> S(3);
  complex<double> a;

  firstIndex  i1;
	secondIndex i2; 
	

	A(2,2) = (9.0,8.0);
	cout << "A(2,2) " << A(2,2) << endl;
	
  real(A(2,3))=1.0;
  imag(A(2,3))=2.0;
  real(A(7,3))=2.0;
  imag(A(7,3))=15.0;
//  real(A(5,0))=30.0;
 // real(A(6,0))=9.0;

  real(B(2,3))=4.0;
  imag(B(2,3))=3.0;
  real(B(7,3))=2.0;
  imag(B(7,3))=1.0;
  
  A(0,0) = complex<double>(1.0,9.0);
  complex<double> x,y,z;
  
 
	double q = 2.0;
  B = complex<double>(q,0)*(1.0*i1)* A;
    // cout << A << B << endl;
	
	cout << "my_pow " << my_pow(3.0,0) << " " << my_pow(3.0,2) << " " << my_pow(3.0,4) << endl;
	cout << "pow(3.0,0) " << pow(3.0,0) << endl;

	
//  global_file << sum(A*conj(B)) << real(sum(A*conj(B))) << endl; 
  /*
  int i=-3; int j=4;
//  cout << C << endl;
  
  Array<complex<double>,1> *S;
  
	S = new Array<complex<double>,1>(3);
  
  real((*S)(0)) = 2.0; imag((*S)(0)) = 3.0;
//  cout << S << endl;
  
//  cout << (*S) << endl;

	delete S;
	delete &C;
	// cout << (*S) << endl;
	*/
} 




