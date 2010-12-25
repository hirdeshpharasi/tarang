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
/*! \file fourier.cc 
 * 
 * @sa fourier.h
 * 
 * @author  M. K. Verma
 * @version 4.0 Parallel 
 * @date	August 2008
 * @bug		none known
 */ 


#include "four_tr.h"

using namespace blitz ;


//*********************************************************************************************

void Init_fftw_plan_FOUR(int N[], Array<complx,3> A)
{
	
	Array<complx,1> temp_row(N[2]);
 	Array<complx,2> temp_plane(N[1],(N[3]/2)+1);
	
	if (N[2] > 1)
	{
		r2c_plan_FOUR = fftw_mpi_plan_dft_r2c_3d(N[1], N[2], N[3], 
							reinterpret_cast<double*>(A.data()), 
							reinterpret_cast<fftw_complex*>(A.data()),  
							MPI_COMM_WORLD, FFTW_MEASURE);
								
		c2r_plan_FOUR = fftw_mpi_plan_dft_c2r_3d(N[1], N[2], N[3], 
							reinterpret_cast<fftw_complex*>(A.data()), 
							reinterpret_cast<double*>(A.data()), 
							MPI_COMM_WORLD, FFTW_MEASURE);	
		
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
	}
	
	else if (N[2] ==1)
	{
		static Array<complx,2> A2d(N[1], (N[3]/2)+1);
		
		r2c_plan_FOUR = fftw_mpi_plan_dft_r2c_2d(N[1], N[3], 
							reinterpret_cast<double*>(A2d.data()), 
							reinterpret_cast<fftw_complex*>(A2d.data()),  
							MPI_COMM_WORLD, FFTW_MEASURE);
		
		c2r_plan_FOUR = fftw_mpi_plan_dft_c2r_2d(N[1], N[3], 
							reinterpret_cast<fftw_complex*>(A2d.data()), 
							reinterpret_cast<double*>(A2d.data()), 
							MPI_COMM_WORLD, FFTW_MEASURE);
	}
	cout << "in plan here 6 " << endl;
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


void ArrayFFTW_FOUR_transpose_order
(
 int N[], 
 Array<complx,3> Atr, 
 Array<complx,3> A
)
{
	if (N[2] > 1)
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

}




//*********************************************************************************************

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

//*********************************************************************************************
//*********************************************************************************************
//*********************************************************************************************
//*********************************************************************************************

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

//*********************************************************************************************

void Zero_pad_last_plane_FOUR(int N[],  Array<complx,3> A)
{
		A(Range(0,local_N1-1),Range(0,N[2]-1),N[3]/2) = 0.0;
}

//*********************************************************************************************

void Norm_FOUR(int N[], Array<complx,3> A) 
{
	A = A/(DP(N[1]) *  DP(N[2]) * DP(N[3]));
}


//*********************************************************************************************

void Xderiv_FOUR(int N[],Array<complx,3> A, Array<complx,3> B, DP kfactor[])
{
	int k1;
	
	for (int l1 = 0; l1 < local_N1; l1++) 			// l1 is the local array-index along x
	{
		k1 = Get_kx_FOUR(l1, N);
		
		B(l1,Range::all(),Range::all()) = 
		    complex<DP>(0, kfactor[1])*((DP) 1.0*k1)*( A(l1,Range::all(),Range::all()) ); 	
	}
}	


//*********************************************************************************************


void Yderiv_FOUR(int N[], Array<complx,3> A, Array<complx,3> B, DP kfactor[])
{
	secondIndex l2;
	B(Range::all(),Range(0,N[2]/2),Range::all()) = 
		complex<DP>(0, kfactor[2])*((DP) 1.0*l2)
		* ( A(Range::all(),Range(0,N[2]/2),Range::all()) ); 
		// l2 = 0:N2/2; k2 = l2.
	
	if (N[2] > 1)
		B(Range::all(),Range(N[2]/2+1,N[2]-1),Range::all()) = 
			complex<DP>(0, kfactor[2])*((DP) 1.0*(l2+1-N[2]/2))*
			( A(Range::all(),Range(N[2]/2+1,N[2]-1),Range::all()) );
			// l2 = 0:N2/2;  k2 = l2-N2/2+1.
}

//*********************************************************************************************


void Zderiv_FOUR(int N[],Array<complx,3> A, Array<complx,3> B, DP kfactor[])
{
	thirdIndex l3;
	B = complex<DP>(0, kfactor[3])*((DP) 1.0*l3)*(A);
}

//********************************	End of four_tr.cc *****************************************


