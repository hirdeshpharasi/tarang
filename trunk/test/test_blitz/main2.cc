

#include <iostream> 
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <complex>
#include <cmath>

BZ_USING_NAMESPACE(blitz)

int main()
{
	int i;
	
	Array<double,1>  A(5);
	
	A = 0;
	
	cout << A;
	cout << " pow  " <<  pow(0.0,1)*5.0 << " " << pow(0.0,0)*5.0 << endl;
	
	 Array<complex<double>,2> B(8,5);
	B(2,2) = (9.0,8.0);
	real(B(2,2)) = 9.0;  imag(B(2,2)) = 8.0;
	cout << "B(2,2) " << B(2,2) << endl;
 
} 





