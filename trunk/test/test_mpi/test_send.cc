
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


#include <iostream> 
#include <mpi.h> 

// using namespace mpi;
using namespace std;

int main(int argc, char *argv[]) 
{ 
	int my_id, numprocs;
	const int master_id = 0;
	const int N = 512;
	int tag1 = 111;
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	
	Array<double,3> A(N, N, N);
	int data_size = N * N * N/numprocs;
	int Nreduced = N/numprocs;
	
	Array<double,3> B(Nreduced, N, N);
	
	A = 1.0;
	B=1.0;
	
	if (my_id == master_id) 
	{
		
		for (int dest = 1; dest <= numprocs-1; dest++) 
		{
			MPI_Send( reinterpret_cast<double*>((B).data()), data_size, 
						MPI_DOUBLE, dest, tag1, MPI_COMM_WORLD );
			cout << "sent " << dest << endl;																										
		}
	}
	
	else		// process is not master
	{
		MPI_Recv( reinterpret_cast<double*>(B.data()), data_size, 
						MPI_DOUBLE, master_id, tag1, MPI_COMM_WORLD, &status);	
						
		cout << "recd " << my_id << endl;				
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Finalize();	
} 




