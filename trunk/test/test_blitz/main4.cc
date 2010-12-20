
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <random/uniform.h> 
#include <random/normal.h> 
#include <time.h> 

#ifdef BZ_HAVE_STD
#include <fstream>
#else
#include <fstream.h>
#endif

BZ_USING_NAMESPACE(blitz)

using namespace blitz;

inline double my_pow(double x, int n) {		// n=0,2,4 only
	if (n==0) return 1.0;
	else if (n==2) return (x*x);
	else if (n==4) return (x*x*x*x);
}

int main()
{
  	cout << "HI "<< endl;
	cout << "my_pow " << my_pow(3.0,0) << " " << my_pow(3.0,2) << " " << my_pow(3.0,4) << endl;
	cout << "pow(3.0,0) " << pow(3.0,0) << endl;
} 