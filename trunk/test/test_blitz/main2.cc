

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
 
} 





