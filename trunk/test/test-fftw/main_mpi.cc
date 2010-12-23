

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
 int local_n0, local_0_start;
  
 
fftw_plan r2c_plan_FOUR, c2r_plan_FOUR;

DP Local_abs_sqr_FOUR(int N[], Array<complx,3> A);
DP Total_abs_sqr_FOUR(int N[], Array<complx,3> A);
	
int main(int argc, char **argv)
{
         const int N0 =4, N1 = 4, N2 = 4;
		 int NN[4]; 
         fftw_plan plan, c2c_2d_forward_plan_FOUR;
	//	 Array<complex<double>,3>  *A;
      //   fftw_complex *data;
         int alloc_local, i, j;
		 int local_n1, local_1_start;
		 DP tot;
		 
		 Array<complex<double>,3>  *A;
    Array<complex<double>,3>  *B;
	 
		NN[1] = N0; NN[2] = N1; NN[3] = N2;
         MPI_Init(&argc, &argv);
		 MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
		 MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		 
		 cout << "MPI " << my_id << " "  << numprocs << endl;
	
         fftw_mpi_init();
     
         /* get local data size and allocate */
         alloc_local = fftw_mpi_local_size_3d(N0, N1, N2, MPI_COMM_WORLD,
                                              &local_n0, &local_0_start);
   //      data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
     
		cout << "Myid , Localn0, localnostart:   " << my_id << " " << local_n0 <<" "  << local_0_start << endl;
		
		A = new Array<complex<double>,3>(local_n0, N1, N2/2+1);
		
		B = new Array<complex<double>,3>(local_n0, N1, N2/2+1);
		
	
	  r2c_plan_FOUR = fftw_mpi_plan_dft_r2c_3d(NN[1], NN[2], NN[3], reinterpret_cast<double*>((*A).data()), 
			     reinterpret_cast<fftw_complex*>((*A).data()),  MPI_COMM_WORLD, FFTW_MEASURE);
	   c2r_plan_FOUR = fftw_mpi_plan_dft_c2r_3d(NN[1], NN[2], NN[3], reinterpret_cast<fftw_complex*>((*A).data()), 
			     reinterpret_cast<double*>((*A).data()), MPI_COMM_WORLD, FFTW_MEASURE);

	c2c_2d_forward_plan_FOUR = fftw_plan_dft_2d(NN[2], NN[3],
												reinterpret_cast<fftw_complex*>((*A).data()), 
												reinterpret_cast<fftw_complex*>((*A).data()),	
												FFTW_FORWARD, FFTW_MEASURE);
	
	  if (my_id == 0) 
			real((*A)(0,0,2)) = 1.0; 
		
		tot = sum(sqr(abs(*A)));
			cout << "Total Before IFFT: myid, tot  " <<  my_id << " " << tot << endl << endl << endl;
	
			 cout << "myid, A   " << my_id << " " << *A << endl;
		//		cout << "myid, B " << my_id << " " << *B << endl;
			
		cout << "Start IFFT " << endl;
       
		
        fftw_execute_dft_c2r(c2r_plan_FOUR, reinterpret_cast<fftw_complex*>((*A).data()), 
		       reinterpret_cast<double*>((*A).data()));
			   
			tot = sum(sqr(abs(*A)));
		cout << " AFTER IFFT: myid,tot " <<  my_id << " " << tot << endl << endl << endl;	   
		 
		   fftw_execute_dft_r2c(r2c_plan_FOUR, reinterpret_cast<double*>((*A).data()), 
			     reinterpret_cast<fftw_complex*>((*A).data()));
     
		
		int data_size = local_n0 * NN[2]* NN[3];
		
		
		 cout << "AFTER FFT, myid, A " << my_id <<  *A << endl;
		 
		 
		 double C[10], D[10];
		
		 for (int i=0; i< 10; i++)  C[i] = 0.0;
		  
		  C[5] = 4.0;
		  
		Array<double, 1> E(10), F(10);
		E=0.0; F=0.0;
		E(1) =  4.0;
		  
		// data_size = local_n0 * NN[2] *NN[3]; 
		MPI_Reduce(E.data(), F.data(), 10, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 

		if (my_id == 0) {
		
			cout << "F "  << F  << endl;
			// for (int i=0; i< 10; i++) 
			//	cout << i << " "  << C[i] << " "  << D[i] << endl;
			
		}
		
		double tot2 = Total_abs_sqr_FOUR(NN, *A);
		if (my_id == 0) 
		 cout << "myid, tot   " << my_id << "   " << tot2 << endl;
			
		 fftw_destroy_plan(r2c_plan_FOUR);
		 fftw_destroy_plan(c2r_plan_FOUR);	
         MPI_Finalize();
}


DP Local_abs_sqr_FOUR(int N[], Array<complx,3> A)
{
//	DP tot = sum(sqr(abs(A))) - sum(sqr(abs(A(Range::all(),Range::all(),0))))/2;

	DP tot = sum(sqr(abs(A)));
	
//	if (local_0_start == 0) 
//		tot = tot -  pow2(abs(A(0,0,0)))/2;
	
	cout << "procs " << my_id << " "  << tot << endl;
	return tot;
                       
  // Subtract |A(k)|^2 of the part of Nz=0 for double counting
  // and |A(0,0)|^2 for mean value
}

DP Total_abs_sqr_FOUR(int N[], Array<complx,3> A)
{
	DP local_total;
	local_total = Local_abs_sqr_FOUR(N, A);
	
	DP total;
	MPI_Reduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, master_id, MPI_COMM_WORLD);	

	if (my_id == master_id) 
		return total;
}



