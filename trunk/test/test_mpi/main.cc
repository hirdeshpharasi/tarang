
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


#include <iostream.h> 
#include <stdio.h> 
#include "mpi.h" 

 
int main(int argc, char *argv[]) 
{ 
    int bufsize, *buf, count; 
    char filename[128]; 
    MPI::Status status; 
 
    MPI::Init(); 
    int myrank = MPI::COMM_WORLD.Get_rank(); 
    int numprocs = MPI::COMM_WORLD.Get_size(); 
    MPI::File thefile = MPI::File::Open(MPI::COMM_WORLD, "testfile", 
                                        MPI::MODE_RDONLY, 
                                        MPI::INFO_NULL); 
    MPI::Offset filesize = thefile.Get_size();  // in bytes 
    filesize    = filesize / sizeof(int);    // in number of ints 
    bufsize     = filesize / numprocs + 1;   // local number to read 
    buf = (int *) malloc (bufsize * sizeof(int)); 
    thefile.Set_view(myrank * bufsize * sizeof(int), 
		     MPI_INT, MPI_INT, "native", MPI::INFO_NULL); 
	
    thefile.Read(buf, bufsize, MPI_INT, &status); 
    count = status.Get_count(MPI_INT); 
    cout << "process " << myrank << " read " << count << " ints" 
	 << endl; 
    thefile.Close(); 
    MPI::Finalize(); 
    return 0; 
} 




