/* Tarang-4.0
 *
 * Copyright (C) 2008, 2009  Mahendra K. Verma
 *
 * Mahendra K. Verma
 * Indian Institute of Technology, Kanpur-208016
 * UP, India
 *
 * mkv@iitk.ac.in
 *
 * This file is part of Tarang-4.0 .
 *
 * Tarang-4.0 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * Tarang-4.0 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Tarang-4.0; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */

/*! \file scft_tr.cc 
 * 
 * @sa scft_tr.h
 * 
 * @author  M. K. Verma
 * @version 4.0
 * @date	August 2008
 * @bug		No known bugs
 */ 

#include "scft_tr.h"

// using namespace blitz ;


/**********************************************************************************************

		Creates fftw_plans given above
		and  put them in the global variables

***********************************************************************************************/

void Init_fftw_plan_SCFT(int NN[], Array<complx,3> A)
{
	Array<DP,1> temp_row(NN[1]);							// temp row for SIN/COS transform
	Array<complx,2> temp_plane(NN[2],(NN[3]/2)+1);			// temp Plane for fourier transform
	Array<complx,1> temp_column_2d((NN[3]/2)+1);
	
	sintr_plan_SCFT = fftw_plan_r2r_1d(NN[1],(DP*)(temp_row.data()), (DP*)(temp_row.data()), 
							FFTW_RODFT10, FFTW_MEASURE);
	costr_plan_SCFT = fftw_plan_r2r_1d(NN[1],(DP*)(temp_row.data()), (DP*)(temp_row.data()), 
							FFTW_REDFT10, FFTW_MEASURE); 
  
	isintr_plan_SCFT = fftw_plan_r2r_1d(NN[1],(DP*)(temp_row.data()), (DP*)(temp_row.data()), 
							FFTW_RODFT01, FFTW_MEASURE);
	icostr_plan_SCFT = fftw_plan_r2r_1d(NN[1],(DP*)(temp_row.data()), (DP*)(temp_row.data()), 
							FFTW_REDFT01, FFTW_MEASURE);  
  
	r2c_plan_SCFT = fftw_plan_dft_r2c_2d(NN[2],NN[3],(DP*)(temp_plane.data()), 
                             (fftw_complex*)(temp_plane.data()),  FFTW_MEASURE);
	c2r_plan_SCFT = fftw_plan_dft_c2r_2d(NN[2],NN[3],(fftw_complex*)(temp_plane.data()), 
                             (DP*)(temp_plane.data()),FFTW_MEASURE);  
	
	r2c_1d_plan_SCFT = fftw_plan_dft_r2c_1d(NN[3], (DP*)(temp_column_2d.data()), 
											(fftw_complex*)(temp_column_2d.data()), FFTW_MEASURE);
	c2r_1d_plan_SCFT = fftw_plan_dft_c2r_1d(NN[3],(fftw_complex*)(temp_column_2d.data()), 
											(DP*)(temp_column_2d.data()), FFTW_MEASURE); 
}


//*********************************************************************************************


void Zero_pad_lastplane_SCFT(int N[], Array<complx,3> A)
{
	A(Range(0,local_N1-1),Range(0,N[2]-1),N[3]/2) = 0.0;
}



///*********************************************************************************************

void Norm_SCFT(int N[], Array<complx,3> A) 
{  
	A = A/(DP(2*N[1]*N[2]*N[3])); 
}


//*********************************************************************************************

void Sintr_row_SCFT(fftw_plan sintr_plan_SCFT, int N[], Array<DP,1> Row)
{
	fftw_execute_r2r(sintr_plan_SCFT, reinterpret_cast<DP*>(Row.data()),
			reinterpret_cast<DP*>(Row.data()));  
}

void Costr_row_SCFT(fftw_plan costr_plan_SCFT, int N[], Array<DP,1> Row)
{
	fftw_execute_r2r(costr_plan_SCFT, reinterpret_cast<DP*>(Row.data()),
			 reinterpret_cast<DP*>(Row.data())); 
}

void ISintr_row_SCFT(fftw_plan isintr_plan_SCFT, int N[], Array<DP,1> Row)
{
	fftw_execute_r2r(isintr_plan_SCFT, reinterpret_cast<DP*>(Row.data()),
			 reinterpret_cast<DP*>(Row.data()));  
}			 

void ICostr_row_SCFT(fftw_plan icostr_plan_SCFT, int N[], Array<DP,1> Row)
{			 
	fftw_execute_r2r(icostr_plan_SCFT, reinterpret_cast<DP*>(Row.data()),
			reinterpret_cast<DP*>(Row.data()));  		 
}	


void FT_Plane_SCFT(fftw_plan r2c_plan_SCFT, int N[], Array<complx,2> Plane)	
{						
	fftw_execute_dft_r2c(r2c_plan_SCFT, reinterpret_cast<DP*>(Plane.data()), 
                           reinterpret_cast<fftw_complex*>(Plane.data()));	
}	

	
void IFT_Plane_SCFT(fftw_plan c2r_plan_SCFT, int N[], Array<complx,2> Plane)
{				
	fftw_execute_dft_c2r(c2r_plan_SCFT, reinterpret_cast<fftw_complex*>(Plane.data()), 
                         reinterpret_cast<DP*>(Plane.data()));
}	

void FT_column_1d_SCFT(fftw_plan r2c_1d_plan_SCFT, int N[], Array<complx,1> temp_column_2d)	
{						
	fftw_execute_dft_r2c(r2c_1d_plan_SCFT, reinterpret_cast<DP*>(temp_column_2d.data()), 
						 reinterpret_cast<fftw_complex*>(temp_column_2d.data()));	
}	


void IFT_column_1d_SCFT(fftw_plan c2r_1d_plan_SCFT, int N[], Array<complx,1> temp_column_2d)
{				
	fftw_execute_dft_c2r(c2r_1d_plan_SCFT, 
						 reinterpret_cast<fftw_complex*>(temp_column_2d.data()), 
                         reinterpret_cast<DP*>(temp_column_2d.data()));
}
							 	 

//*********************************************************************************************


void ArrayShiftRight_SCFT(int N[], Array<complx,3> A)
{
	static Array<complx,2> temp_plane(N[2], (N[3]/2)+1);
	int data_size = 2* N[2] * (N[3]/2+1);
	int tag = 222;
	
	// All the procs except the last one send the last column  to the next-right processor
	if (my_id <= numprocs-2) 
	{												
		temp_plane = A(local_N1-1, Range::all(), Range::all());
		MPI_Send( reinterpret_cast<double*>(temp_plane.data()), data_size, MPI_DOUBLE, 
						my_id+1, tag, MPI_COMM_WORLD );	
	}
						
	for (int i = local_N1-2;  i >= 0; i--)   
	{ 
		temp_plane = A(i,Range::all(),Range::all());
		A(i+1,Range::all(),Range::all()) = temp_plane;
	}
    
	
	if (my_id > 0) 
	{
		MPI_Recv( reinterpret_cast<double*>(temp_plane.data()), data_size, 
						MPI_DOUBLE,my_id-1, tag, MPI_COMM_WORLD, &status );							
		A(0,Range::all(),Range::all()) = temp_plane;
	}
	else		// my_id = 0
		A(0,Range::all(),Range::all()) = 0.0;
		
	MPI_Barrier(MPI_COMM_WORLD);				// Sync the procs
}

//
//

void ArrayShiftLeft_SCFT(int N[], Array<complx,3> A)
{
	static Array<complx,2> temp_plane(N[2], (N[3]/2)+1);
	int data_size = 2* N[2] * (N[3]/2+1);
	int tag = 111;
	
	// All the procs except the last one send the last column  to the next-left processor
	if (my_id > 0) 
	{												
		temp_plane = A(0, Range::all(),Range::all());
		MPI_Send( reinterpret_cast<double*>(temp_plane.data()), data_size, MPI_DOUBLE, 
						my_id-1, tag, MPI_COMM_WORLD );	
	}
						
	for (int i = 1; i <= local_N1-1; i++)   
	{ 
		temp_plane = A(i,Range::all(),Range::all());
		A(i-1,Range::all(),Range::all()) = temp_plane;
	}
	
	if (my_id < numprocs-1) 
	{
		MPI_Recv( reinterpret_cast<double*>(temp_plane.data()), data_size, 
						MPI_DOUBLE,my_id+1, tag, MPI_COMM_WORLD, &status );							
		A(local_N1-1,Range::all(),Range::all()) = temp_plane;
	}
	else		// my_id = numprocs-1
		A(local_N1-1,Range::all(),Range::all()) = 0.0;
		
	MPI_Barrier(MPI_COMM_WORLD);			// Sync the procs
}


/**********************************************************************************************

	SFT(Atr) = A; Atr is in transposed order

	Note: z= N[3]/2 plane does contain any real data (zero everywhere).
	
***********************************************************************************************/

void ArraySFT_SCFT_transpose_order
(
	fftw_plan sintr_plan_SCFT, 
	fftw_plan r2c_plan_SCFT, 
	fftw_plan r2c_1d_plan_SCFT,
	int N[], 
	Array<complx,3> Atr, 
	Array<complx,3> A
)
{
	if (N[2] > 1)
	{
		int l1, l3;
  
		static Array<DP,1> temp_col(N[1]);						// Temp location to hold *Vi(:,j,k)
		static Array<complx,2> temp_plane(N[2], (N[3]/2)+1);     // Temp location to hold *Vi(i,:)
			
		// Array Atr is transposed N2, N1, N3
		// Sin transform along x for all y,z of Atr
		for (l1=0; l1<local_N2; l1++)
			for (l3 = 0; l3 < N[3]/2; l3++)   
			{
				temp_col = real(Atr(l1, Range::all(), l3));  
				Sintr_row_SCFT(sintr_plan_SCFT, N, temp_col);  
				real(Atr(l1, Range::all(), l3)) = temp_col;
		
				temp_col = imag(Atr(l1, Range::all(), l3)); 
				Sintr_row_SCFT(sintr_plan_SCFT, N, temp_col);
				imag(Atr(l1, Range::all(), l3)) = temp_col;
			}
	  
		Zero_pad_lastplane_SCFT(N, Atr);							 // Zero_pad the row=N[2]/2
		Inverse_transpose_array(N, Atr, A);
	  
		// Array A is normal ordered
		// FFT along y-z planes of  A
		for (l1 = 0; l1 < local_N1; l1++)   
		{					
			temp_plane = A(l1, Range::all(), Range::all());
			FT_Plane_SCFT(r2c_plan_SCFT, N, temp_plane);
			A(l1, Range::all(), Range::all()) = temp_plane;
		}
	  
		Norm_SCFT(N, A);
		ArrayShiftRight_SCFT(N, A);	
	}
	
	else if (N[2] == 1)
	{
		int l1;
		
		static Array<DP,1> temp_col(N[1]);					// Temp location to hold *Vi(:,j)
		static Array<complx,1> temp_column_2d((N[3]/2)+1);  // Temp location to hold *Vi(i,:)
		
			// Array Atr is transposed N2, N1, N3
			// Sin transform along x for all y,z of Atr
		for (l1=0; l1<local_N2; l1++)
		{
			temp_col = real(Atr(l1, 0, Range::all()));  
			Sintr_row_SCFT(sintr_plan_SCFT, N, temp_col);  
			real(Atr(l1, 0, Range::all())) = temp_col;
			
			temp_col = imag(Atr(l1, 0, Range::all())); 
			Sintr_row_SCFT(sintr_plan_SCFT, N, temp_col);
			imag(Atr(l1, 0, Range::all())) = temp_col;
		}
		
		Zero_pad_lastplane_SCFT(N, Atr);							 // Zero_pad the row=N[2]/2
		Inverse_transpose_array(N, Atr, A);
		
			// Array A is normal ordered
			// FFT along y-z planes of  A
		for (l1 = 0; l1 < local_N1; l1++)   
		{					
			temp_column_2d = A(l1, 0, Range::all());
			FT_column_1d_SCFT(r2c_1d_plan_SCFT, N, temp_column_2d);
			A(l1, 0, Range::all()) = temp_column_2d;
		}
		
		Norm_SCFT(N, A);
		ArrayShiftRight_SCFT(N, A);
	}
}

//
//
void ArraySFT_SCFT
(
	fftw_plan sintr_plan_SCFT, 
	fftw_plan r2c_plan_SCFT, 
	fftw_plan r2c_1d_plan_SCFT,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_r
)
{
	Transpose_array(N, A, temp_r);
	ArraySFT_SCFT_transpose_order(sintr_plan_SCFT, r2c_plan_SCFT, r2c_1d_plan_SCFT, N, temp_r, A);
}

/**********************************************************************************************

	 CFT(Atr) = A; Atr is in transposed order.  Atr unchanged
	 
	 Note: z= N[3]/2 plane does contain any real data (zero everywhere).

***********************************************************************************************/

void ArrayCFT_SCFT_transpose_order
(
	fftw_plan costr_plan_SCFT, 
	fftw_plan r2c_plan_SCFT, 
	fftw_plan r2c_1d_plan_SCFT,
	int N[], 
	Array<complx,3> Atr, 
	Array<complx,3> A
)
{
	if (N[2] > 1)
	{	
		int l1, l3;
	  
		static Array<DP,1> temp_col(N[1]);						// Temp location to hold *Vi(:,j,k)
		static Array<complx,2> temp_plane(N[2], (N[3]/2)+1);     // Temp location to hold *Vi(i,:) 
										
		// Array Atr is transposed N2, N1, N3
		// Sin transform along x for all y,z of Atr
		for (l1=0; l1<local_N2; l1++)
			for (l3 = 0; l3 < N[3]/2; l3++)   
			{
				temp_col=real(Atr(l1, Range::all(), l3));  
				Costr_row_SCFT(costr_plan_SCFT, N, temp_col);  
				real(Atr(l1, Range::all(), l3))=temp_col;
		
				temp_col=imag(Atr(l1, Range::all(), l3)); 
				Costr_row_SCFT(costr_plan_SCFT, N, temp_col);
				imag(Atr(l1, Range::all(), l3))=temp_col;
			}
		
		Zero_pad_lastplane_SCFT(N, Atr);							// Zero_pad the k=N[3]/2
		Inverse_transpose_array(N, Atr, A);
		
		// Array A is normal ordered
		// FFT along y-z planes of  A
		for (l1=0; l1 < local_N1; l1++)   
		{								
			temp_plane = A(l1, Range::all(), Range::all());
			FT_Plane_SCFT(r2c_plan_SCFT, N, temp_plane);
			A(l1, Range::all(), Range::all()) = temp_plane;
		}

		Norm_SCFT(N,A);
	}
	
		// 2D case
	else if (N[2] == 1)
	{
		int l1;
		
		static Array<DP,1> temp_col(N[1]);					// Temp location to hold *Vi(:,k)
		static Array<complx,1> temp_column_2d((N[3]/2)+1);  // Temp location to hold *Vi(i,:) 
		
		for (l1=0; l1<local_N2; l1++)
		{
			temp_col = real(Atr(l1, 0, Range::all())); 
			Costr_row_SCFT(costr_plan_SCFT, N, temp_col);  
			real(Atr(l1, 0, Range::all())) = temp_col;
			
			temp_col = imag(Atr(l1, 0, Range::all())); 
			Costr_row_SCFT(costr_plan_SCFT, N, temp_col);
			imag(Atr(l1, 0, Range::all())) = temp_col;
		}
		
		Zero_pad_lastplane_SCFT(N, Atr);							// Zero_pad the k=N[3]/2
		Inverse_transpose_array(N, Atr, A);
		
			// Array A is normal ordered
			// FFT along z column of  A
		for (l1=0; l1 < local_N1; l1++)   
		{								
			temp_column_2d = A(l1, 0, Range::all());
			FT_column_1d_SCFT(r2c_1d_plan_SCFT, N, temp_column_2d);
			A(l1, 0, Range::all()) = temp_column_2d;
		}
	}
}

//
//
void ArrayCFT_SCFT
(
	fftw_plan costr_plan_SCFT, 
	fftw_plan r2c_plan_SCFT, 
	fftw_plan r2c_1d_plan_SCFT,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_r
)
{
	Transpose_array(N, A, temp_r);
	ArrayCFT_SCFT_transpose_order(costr_plan_SCFT, r2c_plan_SCFT, r2c_1d_plan_SCFT, N, temp_r, A);
}									

/**********************************************************************************************

	Inverse SFT (A)  = Atr 
	IFT along perp dirn and SIN transform along x dirn of A 

***********************************************************************************************/

void ArrayISFT_SCFT_transpose_order
(
	fftw_plan isintr_plan_SCFT, 
	fftw_plan c2r_plan_SCFT, 
	fftw_plan c2r_1d_plan_SCFT,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> Atr
)  
{
   
	if (N[2] > 1)
	{
		int l1, l3;
		static Array<DP,1> temp_col(N[1]);						// Temp location to hold *Vi(:,j,k)
		static Array<complx,2> temp_plane(N[2], (N[3]/2)+1);     // Temp location to hold *Vi(i,:) 
		
		ArrayShiftLeft_SCFT(N, A);

		// Array A is normal  ordered
		// IFFT along y-z planes of  A
		for (l1=0; l1 < local_N1; l1++)   								// FFT along y
		{
			temp_plane = A(l1, Range::all(), Range::all());
			IFT_Plane_SCFT(c2r_plan_SCFT, N, temp_plane);
			A(l1, Range::all(), Range::all()) = temp_plane;
		}
		
		Transpose_array(N, A, Atr);
	  
		// Array Atr is transposed N2, N1, N3
		// ISin transform along x for all y,z of Atr
		for (l1=0; l1 < local_N2; l1++)											// Sin tr alog x
			 for (l3=0; l3<N[3]/2; l3++)  
			 {        
				temp_col = real(Atr(l1, Range::all(), l3)); 
				ISintr_row_SCFT(isintr_plan_SCFT, N, temp_col);  
				real(Atr(l1, Range::all(), l3)) = temp_col;
		
				temp_col = imag(Atr(l1, Range::all(), l3)); 
				ISintr_row_SCFT(isintr_plan_SCFT, N, temp_col);  
				imag(Atr(l1, Range::all(), l3)) = temp_col;
			}
			
		Zero_pad_lastplane_SCFT(N, Atr);	
	}
	
		// 2D case..
	else if (N[2] == 1)
	{
		int l1;
		static Array<DP,1> temp_col(N[1]);						// Temp location to hold *Vi(:,j,k)
		static Array<complx,1> temp_column_2d((N[3]/2)+1);     // Temp location to hold *Vi(i,:) 
		
		ArrayShiftLeft_SCFT(N, A);
		
			// Array A is normal  ordered
			// IFFT along y-z planes of  A
		for (l1=0; l1 < local_N1; l1++)   								// FFT along y
		{
			temp_column_2d = A(l1, 0, Range::all());
			IFT_column_1d_SCFT(c2r_1d_plan_SCFT, N, temp_column_2d);
			A(l1, Range::all()) = temp_column_2d;
		}
		
		Transpose_array(N, A, Atr);
		
			// Array Atr is transposed N2, N1, N3
			// ISin transform along x for all y,z of Atr
		for (l1=0; l1 < local_N2; l1++)											// Sin tr alog x 
		{        
			temp_col = real(Atr(l1, 0, Range::all())); 
			ISintr_row_SCFT(isintr_plan_SCFT, N, temp_col);  
			real(Atr(l1, 0, Range::all())) = temp_col;
			
			temp_col = imag(Atr(l1, 0, Range::all())); 
			ISintr_row_SCFT(isintr_plan_SCFT, N, temp_col);  
			imag(Atr(l1, 0, Range::all())) = temp_col;
		}
		
		Zero_pad_lastplane_SCFT(N, Atr);
	}
}

//
//
void ArrayISFT_SCFT
(
	fftw_plan isintr_plan_SCFT, 
	fftw_plan c2r_plan_SCFT, 
	fftw_plan c2r_1d_plan_SCFT,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_r
)
{
	ArrayISFT_SCFT_transpose_order(isintr_plan_SCFT, c2r_plan_SCFT, c2r_1d_plan_SCFT, N, A, temp_r);
	Inverse_transpose_array(N, temp_r, A);
}

/**********************************************************************************************

				Inverse CFT(A) = Atr
	IFT along y-z dirn and COS transform along x dirn of A

***********************************************************************************************/

void ArrayICFT_SCFT_transpose_order
(
	fftw_plan icostr_plan_SCFT, 
	fftw_plan c2r_plan_SCFT, 
	fftw_plan c2r_1d_plan_SCFT,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> Atr
)  
{
  
	if (N[2] > 1)
	{	
		int l1, l3;
		static Array<DP,1> temp_col(N[1]);						// Temp location to hold *Vi(:,j,k)
		static Array<complx,2> temp_plane(N[2], (N[3]/2)+1);     // Temp location to hold *Vi(i,:)
		
		// Array A is normal ordered
		// IFFT along y-z planes of  A		
		for (l1=0; l1 < local_N1; l1++)   								// FFT along y-z plane
		{
			temp_plane = A(l1, Range::all(), Range::all());
			IFT_Plane_SCFT(c2r_plan_SCFT, N, temp_plane);
			A(l1, Range::all(), Range::all()) = temp_plane;
		}
		
		Transpose_array(N, A, Atr);
		
		for (l1=0; l1 < local_N2; l1++)									// Cos tr along x
			for (l3=0; l3<N[3]/2; l3++)  
			{					
				temp_col = real(Atr(l1, Range::all(), l3));  
				ICostr_row_SCFT(icostr_plan_SCFT, N, temp_col);  
				real(Atr(l1, Range::all(), l3)) = temp_col;
		
				temp_col = imag(Atr(l1, Range::all(), l3)); 
				ICostr_row_SCFT(icostr_plan_SCFT, N, temp_col);  
				imag(Atr(l1, Range::all(), l3)) = temp_col;
			}
			
		Zero_pad_lastplane_SCFT(N, Atr);
	}
	
		// 2D case
	else if (N[2] == 1)
	{
		int l1;
		static Array<DP,1> temp_col(N[1]);					// Temp location to hold *Vi(:,k)
		static Array<complx,1> temp_column_2d((N[3]/2)+1);  // Temp location to hold *Vi(i,:)
		
			// Array A is normal ordered
			// IFFT along y-z planes of  A		
		for (l1=0; l1 < local_N1; l1++)   								// FFT along y-z plane
		{
			temp_column_2d = A(l1, 0, Range::all());
			IFT_column_1d_SCFT(c2r_1d_plan_SCFT, N, temp_column_2d);
			A(l1, 0, Range::all()) = temp_column_2d;
		}
		
		Transpose_array(N, A, Atr);
	
		for (l1=0; l1 < local_N2; l1++)									// Cos tr along x
		{					
			temp_col = real(Atr(l1, 0, Range::all()));  
			ICostr_row_SCFT(icostr_plan_SCFT, N, temp_col);  
			real(Atr(l1, 0, Range::all())) = temp_col;
			
			temp_col = imag(Atr(l1, 0, Range::all())); 
			ICostr_row_SCFT(icostr_plan_SCFT, N, temp_col);  
			imag(Atr(l1, 0, Range::all())) = temp_col;
		}
		
		Zero_pad_lastplane_SCFT(N, Atr);
		
	}
}

//
//
void ArrayICFT_SCFT
(
	fftw_plan icostr_plan_SCFT, 
	fftw_plan c2r_plan_SCFT,
	fftw_plan c2r_1d_plan_SCFT,
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_r
)
{
	ArrayICFT_SCFT_transpose_order(icostr_plan_SCFT, c2r_plan_SCFT, c2r_1d_plan_SCFT, N, A, temp_r);

	Inverse_transpose_array(N, temp_r, A);
}
						


/**********************************************************************************************

	Derivative along x;  Bk_COS[i]=pi*i1*Ak_SIN[i]; Bk_SIN[]=-pi*i1*Ak_COS[i]

***********************************************************************************************/

void  Xderiv_Sin_SCFT(int NN[], Array<complx,3> A, Array<complx,3> B, DP kfactor[])
{	
	int k1;
	
	for (int l1 = 0; l1 < local_N1; l1++) 			// l1 is the local array-index along x
	{
		k1 = Get_kx_SCFT(l1, NN);
		B(l1,Range::all(),Range::all()) = 
		    complex<DP>(kfactor[1], 0)*((DP) 1.0*k1)*( A(l1,Range::all(),Range::all()) ); 	
	}
}


void Xderiv_Cos_SCFT(int NN[], Array<complx,3> A, Array<complx,3> B, DP kfactor[])
{
	int k1;
	
	for (int l1 = 0; l1 < local_N1; l1++) 			// l1 is the local array-index along x
	{
		k1 = Get_kx_SCFT(l1, NN);
		B(l1,Range::all(),Range::all()) = 
		    -complex<DP>(kfactor[1], 0)*((DP) 1.0*k1)*( A(l1,Range::all(),Range::all()) ); 	
	}
}

/**********************************************************************************************

		Derivative along y;  B(k) = i*ky*A(k)

***********************************************************************************************/

// Note: In the first half- ky=i2;
// In the second half- i2=0:Ny/-1; fftw-index=(Ny/2 +1+i2); FT-index=fftw-index-N=(i2+1-Ny/2)
void Yderiv_SCFT(int NN[], Array<complx,3> A, Array<complx,3> B, DP kfactor[])
{
	secondIndex i2;
  
	B(Range::all(),Range(0,NN[2]/2),Range::all()) = 
		complex<DP>(0,kfactor[2])*((DP) 1.0*i2)
					*( A(Range::all(),Range(0,NN[2]/2),Range::all()) ); 
  
	B(Range::all(),Range(NN[2]/2+1,NN[2]-1),Range::all()) = 
           complex<DP>(0,kfactor[2])*((DP) 1.0*(i2+1-NN[2]/2))
					*( A(Range::all(),Range(NN[2]/2+1,NN[2]-1),Range::all()) );
}

/**********************************************************************************************

		Derivative along z; Derivative along z

***********************************************************************************************/


void Zderiv_SCFT(int NN[], Array<complx,3> A, Array<complx,3> B, DP kfactor[])
{
	thirdIndex i3;
	
	B = complex<DP>(0, kfactor[3])*((DP) 1.0*i3) * (A);     
}


//******************************** End of scft_tr.cc  *****************************************






