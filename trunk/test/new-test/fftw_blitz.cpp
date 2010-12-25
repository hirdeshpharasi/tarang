/*  fftw_performance.cpp

    Performance testing of FFTW-MPI transforms.
    Copyright 2010 Mani Chandra <mchandra@iitk.ac.in>
		           Supriyo Paul <supriyo@iitk.ac.in>
····
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
····
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
····
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
    MA 02110-1301, USA.
*/

#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <time.h>
#include <string>
#include <fstream>
// #include "param.h"

#ifndef DP
#define DP double
#endif

#ifndef complx
#define complx  complex<DP>
#endif

using namespace blitz;

int my_id;								// My process id
int numprocs;							// No of processors
const int master_id=0;						// Id of master proc
ptrdiff_t alloc_local;

int globalvar_fftw_original_switch = 0;

ptrdiff_t local_N1, local_N1_start;		// N1 size and start of i1 in the currentproc
ptrdiff_t local_N2, local_N2_start;

// for fftw_original
fftw_plan r2c_plan_FOUR, c2r_plan_FOUR;

// for split fftw
fftw_plan r2c_2d_plan_FOUR, c2r_2d_plan_FOUR;
fftw_plan c2c_1d_forward_plan_FOUR, c2c_1d_inverse_plan_FOUR;

void FTr2c_plane_FOUR(int N[], Array<complx,2> Plane);
void FTc2r_plane_FOUR(int N[], Array<complx,2> Plane);
void FT_col_FOUR(int N[], Array<complx,1> Col);
void IFT_col_FOUR(int N[], Array<complx,1> Col);


void initialize(int N[], Array<complx,3> A);
void Norm_FOUR(int N[], Array<complx,3> A);
void Zero_pad_last_plane_FOUR(int N[],  Array<complx,3> A);


void ArrayFFT_FOUR
(
 int N[], 
 Array<complx,3> A,
 Array<complx,3> temp_r
 );

void ArrayFFT_FOUR_transpose_order
(
 int N[], 
 Array<complx,3> Atr, 
 Array<complx,3> A
 );

void ArrayFFTW_FOUR_transpose_order
(
 int N[], 
 Array<complx,3> Atr, 
 Array<complx,3> A
 );
void ArrayFFTW_FOUR
(
 int N[], 
 Array<complx,3> A,
 Array<complx,3> temp_r
 );


void ArrayIFFT_FOUR_transpose_order
(
 int N[], 
 Array<complx,3> A, 
 Array<complx,3> Atr
);

void ArrayIFFT_FOUR
(
 int N[], 
 Array<complx,3> A,
 Array<complx,3> temp_r
 );

void Transpose_array(int N[], Array<complx,3> A, Array<complx,3> Atr);
void Inverse_transpose_array(int N[], Array<complx,3> Atr, Array<complx,3> A);

string itostring(int value);

int main(int argc, char** argv)
{
    const int N1 =4, N2 = 4, N3 = 4;
	int N[4]; 
    const int dim = 3;
    const int number_of_transforms = 1;
    unsigned int numBytes;
	
	N[1] = N1; N[2] = N2; N[3] = N3;
   
    time_t start, end;
    DP diff_time;
    
    string outfile_name, numstr1, numstr2, numstr3, ext;
    ofstream outfile;


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    fftw_mpi_init();
	
	local_N1 = N[1]/numprocs;			
	local_N2 = N[2]/numprocs;
	local_N1_start = my_id * local_N1;
	local_N2_start = my_id * local_N2;
	
	cout << "my_id , local_Ni, local_Ni_start " << my_id << " " << local_N1 <<  " " << local_N2
	<< " " << local_N1_start << " " << local_N2_start << endl;

	cout << "NUMPROC " << numprocs << endl << endl;
	
	Array<complx,1> temp_row(N[2]);
 	Array<complx,2> temp_plane(N[1],N[3]/2+1);
	
	Array<complex<double>,3>  A(local_N1, N2, N3/2+1);
    Array<complex<double>,3>  temp_r(local_N2, N1, N3/2+1);

//	alloc_local = fftw_mpi_local_size_3d(N1, N2, N3, MPI_COMM_WORLD, 
//							 &local_N1, &local_N1_start, &local_N2, local_N2_start);
//	numBytes = sizeof(fftw_complex)*alloc_local;
//	data_array = (fftw_complex*) fftw_malloc(numBytes);
	
	
	r2c_plan_FOUR = fftw_mpi_plan_dft_r2c_3d(N[1], N[2], N[3], 
							 reinterpret_cast<double*>(A.data()), 
							 reinterpret_cast<fftw_complex*>(A.data()),  
							 MPI_COMM_WORLD, FFTW_MEASURE);
	
	c2r_plan_FOUR = fftw_mpi_plan_dft_c2r_3d(N[1], N[2], N[3], 
							 reinterpret_cast<fftw_complex*>(A.data()), 
							 reinterpret_cast<double*>(A.data()), 
							 MPI_COMM_WORLD, FFTW_MEASURE);	
	
	/*		r2c_2d_forward_plan_FOUR = fftw_plan_dft_r2c_2d(N[1], N[3], 
	 reinterpret_cast<double*>(temp_plane_forward.data()), 
	 reinterpret_cast<fftw_complex*>(temp_plane_forward.data()),  
	 FFTW_MEASURE); */
	
	r2c_2d_plan_FOUR = fftw_plan_dft_r2c_2d(N[1], N[3], 
								reinterpret_cast<double*>(temp_plane.data()),
								reinterpret_cast<fftw_complex*>(temp_plane.data()),
								 FFTW_MEASURE);
	
	c2r_2d_plan_FOUR = fftw_plan_dft_c2r_2d(N[1], N[3], 
								reinterpret_cast<fftw_complex*>(temp_plane.data()),
								reinterpret_cast<double*>(temp_plane.data()),
								FFTW_MEASURE);
	
	c2c_1d_forward_plan_FOUR = fftw_plan_dft_1d(N[2], 
								reinterpret_cast<fftw_complex*>(temp_row.data()), 
								reinterpret_cast<fftw_complex*>(temp_row.data()), 
								FFTW_FORWARD, FFTW_MEASURE);
	
	c2c_1d_inverse_plan_FOUR = fftw_plan_dft_1d(N[2], 
								reinterpret_cast<fftw_complex*>(temp_row.data()), 
								reinterpret_cast<fftw_complex*>(temp_row.data()), 
								FFTW_BACKWARD, FFTW_MEASURE);

	if (my_id==master_id)
		cout << "Allocated " << " array on " << numprocs << " processors" << endl << endl;

	initialize(N, A);

	if (my_id==master_id)
		cout << "Initialization complete" << endl << "Beginning execution..." << endl << endl;


    time(&start);
	
	cout << "orig " << A << endl << endl << endl << endl << endl;
	sleep(1);
	
	ArrayIFFT_FOUR(N, A, temp_r);
	Zero_pad_last_plane_FOUR(N, A);
	
	cout << "after ifft " << my_id << " " << A(Range::all(), Range::all(), Range(0,N[3]/2-1)) 
		<< endl << endl << endl;
	
	sleep(1);
	
	ArrayFFT_FOUR(N, A, temp_r);
	
	cout << "after ifft & fft " << my_id << A << endl << endl << endl;
	
	sleep(1);
	
//	fftw_execute(forward_plan); 
//  fftw_execute(backward_plan); 

    time(&end);

    if (my_id==master_id)
    {
        cout << "Execution complete" << endl << endl;
        cout << "Number of transforms performed = " << number_of_transforms << endl << endl;
    }

//    fftw_destroy_plan(forward_plan);
//   fftw_destroy_plan(backward_plan);
//    fftw_free(data_array);

    MPI_Finalize();

    if (my_id==master_id)
    {   
/*		numstr1 = itostring(size_of_array);
    	numstr2 = itostring(dim);
    	numstr3 = itostring(number_of_transforms);
    	ext = ".txt";
    	outfile_name = "SIZE_" + numstr1 + "_DIM_" + numstr2 + "_TRANSFORMS_" + numstr3 + ext;
    	outfile.open(outfile_name.c_str());  */


    	diff_time = difftime(end, start);
        
		cout << endl << "Progam terminating here." << endl << endl \
        << "Total time elapsed =  " << diff_time << " sec" << endl;
      /*  
	//	outfile << "Size = " << size_of_array << endl;
		outfile << "Dimension = " << dim << endl;
		outfile << "Number of transforms = " << number_of_transforms << endl;
		outfile << "Number of processors = " << numprocs << endl;
		outfile << "Elapsed Time = " << diff_time << " sec" << endl;

		outfile.close(); */
		cout << "Finished. Quitting..." << endl; 
    }
}

string itostring(int value) 
{
    stringstream sstr;
    sstr << value;
    return sstr.str();
}
//*********************************************************************************************

void FTr2c_plane_FOUR(int N[], Array<complx,2> Plane)
{
	fftw_execute_dft_r2c(r2c_2d_plan_FOUR, reinterpret_cast<DP*>(Plane.data()), 
					 reinterpret_cast<fftw_complex*>(Plane.data()));
}

void FTc2r_plane_FOUR(int N[], Array<complx,2> Plane)
{
	fftw_execute_dft_c2r(c2r_2d_plan_FOUR, reinterpret_cast<fftw_complex*>(Plane.data()), 
					 reinterpret_cast<DP*>(Plane.data()));
}

void FT_col_FOUR(int N[], Array<complx,1> Col)
{
	fftw_execute_dft(c2c_1d_forward_plan_FOUR, reinterpret_cast<fftw_complex*>(Col.data()), 
						 reinterpret_cast<fftw_complex*>(Col.data()));
}

void IFT_col_FOUR(int N[], Array<complx,1> Col)
{
	fftw_execute_dft(c2c_1d_inverse_plan_FOUR, reinterpret_cast<fftw_complex*>(Col.data()), 
					 reinterpret_cast<fftw_complex*>(Col.data()));
}

//*********************************************************************************************

void initialize(int N[], Array<complx,3> A)
{
	A=0.0;
	if (my_id == 0)  {
		real(A(1,1,1))= 5.0;imag(A(1,1,1))= 4.0;
	}
		
	  
}

void Norm_FOUR(int N[], Array<complx,3> A) 
{
	A = A/(DP(N[1]*N[2]*N[3]));  
}

void Zero_pad_last_plane_FOUR(int N[],  Array<complx,3> A)
{
	A(Range(0,local_N1-1),Range(0,N[2]-1),N[3]/2) = 0.0;
}

//**********************************************************************************************
void ArrayFFT_FOUR
(
 int N[], 
 Array<complx,3> A,
 Array<complx,3> temp_r
 )
{
	Zero_pad_last_plane_FOUR(N, A); 
	ArrayFFTW_FOUR(N, A, temp_r); 
	Norm_FOUR(N, A); 
}


void ArrayFFT_FOUR_transpose_order
(
 int N[], 
 Array<complx,3> Atr, 
 Array<complx,3> A
 )
{
	Zero_pad_last_plane_FOUR(N, A); 
	ArrayFFTW_FOUR_transpose_order(N, Atr, A); 
	Norm_FOUR(N, A); 
}

//**********************************************************************************************

void ArrayFFTW_FOUR_transpose_order
(
 int N[], 
 Array<complx,3> Atr, 
 Array<complx,3> A
 )
{

	
	static Array<complx,1> temp_col(N[2]);						// Temp location to hold *Vi(:,j,k)
	static Array<complx,2> temp_plane(N[1], (N[3]/2)+1);     // Temp location to hold *Vi(i,:)
	
	for (int l1=0; l1<local_N2; l1++) {
		temp_plane = Atr(l1, Range::all(), Range::all());
		FTr2c_plane_FOUR(N, temp_plane);
		Atr(l1, Range::all(), Range::all()) = temp_plane;
	}
		
	Inverse_transpose_array(N, Atr, A);
	
	for (int l1=0; l1<local_N1; l1++) 
		for (int l3=0; l3<=N[3]/2; l3++) {
			temp_col = A(l1, Range::all(), l3);
			FT_col_FOUR(N, temp_col);
			A(l1, Range::all(), l3) = temp_col;
		}
}	
//************************************************************************************************

void ArrayFFTW_FOUR
(
 int N[], 
 Array<complx,3> A,
 Array<complx,3> temp_r
 ) 
{
	
	if (globalvar_fftw_original_switch == 1) {
		static Array<complx,2> A2d(N[1], (N[3]/2)+1); 
		
		if (N[2] > 1) {	
			fftw_execute_dft_r2c(r2c_plan_FOUR, reinterpret_cast<double*>(A.data()), 
								 reinterpret_cast<fftw_complex*>(A.data()));
		}
		
		else if (N[2] == 1) {
			A2d(Range::all(), Range::all()) = A(Range::all(), 0, Range::all());
			
			fftw_execute_dft_r2c(r2c_plan_FOUR, reinterpret_cast<double*>(A2d.data()), 
								 reinterpret_cast<fftw_complex*>(A2d.data())); 
			
			A(Range::all(), 0, Range::all()) = A2d(Range::all(), Range::all());
		}
	}
	
	else {
		Transpose_array(N, A, temp_r);
		ArrayFFTW_FOUR_transpose_order(N, temp_r, A);
	}
	
}

//************************************************************************************************

void ArrayIFFT_FOUR_transpose_order
(
 int N[], 
 Array<complx,3> A, 
 Array<complx,3> Atr
 )
{
	if (N[2] > 1)
	{
		
		static Array<complx,1> temp_col(N[2]);						// Temp location to hold *Vi(:,j,k)
		static Array<complx,2> temp_plane(N[1], (N[3]/2)+1);     // Temp location to hold *Vi(i,:)
		

		for (int l1=0; l1<local_N1; l1++) 
			for (int l3=0; l3<=N[3]/2; l3++) {
				temp_col = A(l1, Range::all(), l3);
				IFT_col_FOUR(N, temp_col);
				A(l1, Range::all(), l3) = temp_col;
			}

		Transpose_array(N, A, Atr);

		for (int l1=0; l1 < local_N2; l1++)	{
			temp_plane = Atr(l1, Range::all(), Range::all());
			FTc2r_plane_FOUR(N, temp_plane);
			Atr(l1, Range::all(), Range::all()) = temp_plane;
		}

		Zero_pad_last_plane_FOUR(N, Atr);
	}
}
//***********************************************************************************************

void ArrayIFFT_FOUR
(
 int N[], 
 Array<complx,3> A,
 Array<complx,3> temp_r
 )
{
	
	if (globalvar_fftw_original_switch == 1) {
		static Array<complx,2> A2d(N[1], (N[3]/2)+1);
		
		if (N[2] > 1)
			fftw_execute_dft_c2r(c2r_plan_FOUR, reinterpret_cast<fftw_complex*>(A.data()), 
								 reinterpret_cast<double*>(A.data()));
		
		
		else if (N[2] == 1) {
			A2d(Range::all(), Range::all()) = A(Range::all(), 0, Range::all());
			
			
			fftw_execute_dft_c2r(c2r_plan_FOUR, reinterpret_cast<fftw_complex*>(A2d.data()), 
								 reinterpret_cast<double*>(A2d.data()));
			
			A(Range::all(), 0, Range::all()) = A2d(Range::all(), Range::all());
		}
	}
	
	else {
		ArrayIFFT_FOUR_transpose_order(N, A, temp_r);
		Inverse_transpose_array(N, temp_r, A);
	}
	
}

//***********************************************************************************************

void Transpose_array(int N[], Array<complx,3> A, Array<complx,3> Atr)
{
	
	
	static Array<complx,2> Axy(local_N1, N[2]);	
	static Array<complx,2> Axy_tr(N[2], local_N1);
	
	static Array<complx,2> temp1(N[1], local_N2);
	static Array<complx,2> temp2(N[2], local_N1);
	
	static Array<complx,2> buffer1_piece(local_N1, local_N2);		
	static Array<complx,2> buffer2_piece(local_N2, local_N1);
	
	int data_size = 2* local_N1 * local_N2;								
	// 2 for complex to double

	if (numprocs == 1)
		Atr = A.transpose(1,0,2);
	
	else
		for (int i3=0; i3<=N[3]/2; i3++)
		{
			
			Axy = A(Range::all(), Range::all(), i3);
			Axy_tr = Axy.transpose(1,0);
			
			// The jth block sent from process i is received by process j and 
			// is placed in the ith block of recvbuf. 
			
			MPI_Alltoall(reinterpret_cast<double*>(Axy_tr.data()), data_size, MPI_DOUBLE, 
						 reinterpret_cast<double*>(temp2.data()),
						 data_size,  MPI_DOUBLE, MPI_COMM_WORLD);	
			
			for (int source = 0; source < numprocs; source++) 
			{																	
				buffer2_piece = temp2(Range(source*local_N2,(source+1)*local_N2-1), 
									  Range::all());
				
				buffer1_piece = buffer2_piece.transpose(1,0);				
				
				temp1(Range(source*local_N1,(source+1)*local_N1-1), Range::all()) 
				= buffer1_piece;
			}
			
			Atr(Range::all(), Range::all(), i3) = temp1.transpose(1,0);
		}

}	
	

//**********************************************************************************************

void Inverse_transpose_array(int N[], Array<complx,3> Atr, Array<complx,3> A)
{
	
		
	static Array<complx,2> Axy(N[1], local_N2);
	static Array<complx,2> Atr_xy(local_N2, N[1]);	
	
	static Array<complx,2> temp1(N[1], local_N2);
	static Array<complx,2> temp2(N[2], local_N1);
	
	static Array<complx,2> buffer1_piece(local_N1, local_N2);		
	static Array<complx,2> buffer2_piece(local_N2, local_N1);
	
	int data_size = 2* local_N1 * local_N2;								
	// 2 for complex to double
	
	if (numprocs == 1)
		A = Atr.transpose(1,0,2);
	
	else
		for (int i3=0; i3<=N[3]/2; i3++)
		{
			
			Atr_xy = Atr(Range::all(), Range::all(), i3);
			Axy = Atr_xy.transpose(1,0);
			
			// The jth block sent from process i is received by process j and 
			// is placed in the ith block of recvbuf. 
			
			MPI_Alltoall(reinterpret_cast<double*>(Axy.data()), data_size, MPI_DOUBLE, 
						 reinterpret_cast<double*>(temp1.data()),
						 data_size,  MPI_DOUBLE, MPI_COMM_WORLD);	
			
			for (int source = 0; source < numprocs; source++) 
			{																	
				buffer1_piece = temp1(Range(source*local_N1,(source+1)*local_N1-1), 
									  Range::all());
				
				buffer2_piece = buffer1_piece.transpose(1,0);				
				
				temp2(Range(source*local_N2,(source+1)*local_N2-1), Range::all()) 
				= buffer2_piece;
			}
			
			A(Range::all(), Range::all(), i3) = temp2.transpose(1,0);
		}	

}
	
