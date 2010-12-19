

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

using namespace ranlib; 

inline double my_arg(complex<double> z)
{
	double temp = arg(z);
	return (temp >= 0) ? temp : temp + 2*M_PI; 
}

inline double my_pow(double x, int n) {		// n=0,2,4 only
	if (n==0) return 1.0;
	else if (n==2) return (x*x);
	else if (n==4) return (x*x*x*x);
}


void test_fn();

Uniform<double> x; 

int main() 
{ 
// At start of program, seed with the system time so we get 
// a different stream of random numbers each run. 
//	Uniform<float> x; 
	x.seed((unsigned int)time(0)); 

	double a,b,angle;
	complex<double> y, z, I;
	
	I = complex<double>(0,1);
	for (int i=0; i<=10; i++) {
		
		angle = 2*M_PI*(x.random()-0.5);
		z = exp(I*angle);
		cout <<  z << " " << angle*180/M_PI << endl;
		cout << arg(z) << " "  <<  my_arg(z)*180/M_PI << "  " << angle << endl << endl;
	}	
 
	cout << "HI "<< endl;
	cout << "my_pow " << my_pow(3.0,0) << " " << my_pow(3.0,2) << " " << my_pow(3.0,4) << endl;
	cout << "pow(3.0,0) " << pow(0.0,0) << " " << pow(0.0,2) << endl;
	
	
} 

void test_fn()
{

	float yf;
	
	for (int i=0; i<10; i++) 
	
	cout << "Inside fn " <<  x.random() << endl;
	
}
