

#include <blitz/array.h>
#include <complex>
#include <cmath>

#ifdef SEEK_SET
#undef SEEK_SET
#endif
#ifdef SEEK_CUR
#undef SEEK_CUR
#endif
#ifdef SEEK_END
#undef SEEK_END
#endif

#include <mpi.h>
#include <fftw3-mpi.h>

#ifndef DP
#define DP double
#endif

#ifndef complx
#define complx  complex<DP>
#endif

// #include <iostream>
// #include <blitz/array.h>
// #include <complex>
// #include <cmath>
// #include <fftw3_mpi.h>
//#include "fourier.h"

using namespace blitz;

int my_id;								// My process id
int numprocs;							// No of processors
const int master_id=0;						// Id of master proc  
 	
int main(int argc, char **argv)
{
	
	int root  = 0;
	
	string prog_kind;
	string data_dir;

	const char    *msgbuf1;
	
	int datasize1, datasize2;
	const char *string_array1;
	const char *string_array2;
	ostringstream str_buffer1, str_buffer2;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		 
	
	if (my_id == master_id) 
	{
		prog_kind = "INC_FLUID";  
		data_dir = "HISSS";
		msgbuf1 = "TESTMKV";
		string_array1 = "INC_FLUID";			datasize1 = prog_kind.size();
		string_array2 = "HISSS";				datasize2 = data_dir.size();
		
		// strcpy(prog_kind, string_array1.c_tr());
	}
	

	MPI_Bcast( &datasize1, 1, MPI_INT, root, MPI_COMM_WORLD);
	MPI_Bcast( &datasize2, 1, MPI_INT, root, MPI_COMM_WORLD);
	
//	MPI_Bcast( &msgbuf1, 7, MPI_CHAR, root, MPI_COMM_WORLD);

	
	MPI_Bcast( &string_array1, datasize1, MPI_CHAR, root, MPI_COMM_WORLD);
	MPI_Bcast( &string_array2, datasize2, MPI_CHAR, root, MPI_COMM_WORLD);
	MPI_Bcast( &msgbuf1, 7, MPI_CHAR, root, MPI_COMM_WORLD);





	if (my_id != master_id)
	{
		str_buffer1 << string_array1;	prog_kind = str_buffer1.str(); 
	//	str_buffer.str("");
	//	str_buffer2 << string_array2;	data_dir = str_buffer2.str();
	}

		if (my_id != master_id)
	{
	//	str_buffer1 << string_array1;	prog_kind = str_buffer1.str(); 
	//	str_buffer.str("");
		cout << "HERE 1 " << endl;
		str_buffer1 << string_array2;	data_dir = str_buffer1.str();
	}


	
	cout << "IN Bcast myid " << my_id << " "  << datasize1 << " "  << datasize2 << " " 
		<< msgbuf1 <<  " " <<  string_array1 << " " << string_array2 << endl;
		
//	cout << "IN Bcast2  myid " << my_id << " " << prog_kind << " "  << data_dir << endl;
//	cout << "IN Bcast myid " << my_id << " "  << datasize1 << endl;

//	cout << "HERE 1 " << endl;
	
	 MPI_Finalize();
	 
//	 cout << "HERE 2 " << endl;
}


