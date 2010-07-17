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
}


//*********************************************************************************************


void ArrayFFTW_FOUR(fftw_plan r2c_plan_FOUR, int N[], Array<complx,3> A) 
{
	
	static Array<complx,2> A2d(N[1], (N[3]/2)+1);    
	
	if (N[2] > 1)
	{	
		fftw_execute_dft_r2c(r2c_plan_FOUR, reinterpret_cast<double*>(A.data()), 
			     reinterpret_cast<fftw_complex*>(A.data()));
	}
	
	else if (N[2] == 1)
	{
		
		A2d(Range::all(), Range::all()) = A(Range::all(), 0, Range::all());
		
		fftw_execute_dft_r2c(r2c_plan_FOUR, reinterpret_cast<double*>(A2d.data()), 
			 reinterpret_cast<fftw_complex*>(A2d.data())); 
		
		A(Range::all(), 0, Range::all()) = A2d(Range::all(), Range::all());
	}	
}

//*********************************************************************************************

void ArrayFFT_FOUR(fftw_plan r2c_plan_FOUR, int N[], Array<complx,3> A)
{
	Zero_pad_last_plane_FOUR(N, A); 
	ArrayFFTW_FOUR(r2c_plan_FOUR, N, A); 
	Norm_FOUR(N, A); 
}

//*********************************************************************************************

void ArrayIFFT_FOUR(fftw_plan c2r_plan_FOUR, int N[], Array<complx,3> A)
{
	static Array<complx,2> A2d(N[1], (N[3]/2)+1); 
	
	if (N[2] > 1)
	{	
		fftw_execute_dft_c2r(c2r_plan_FOUR, reinterpret_cast<fftw_complex*>(A.data()), 
			reinterpret_cast<double*>(A.data()));
	}
	
	else if (N[2] == 1)
	{
		A2d(Range::all(), Range::all()) = A(Range::all(), 0, Range::all());
	
		
		fftw_execute_dft_c2r(c2r_plan_FOUR, reinterpret_cast<fftw_complex*>(A2d.data()), 
							 reinterpret_cast<double*>(A2d.data()));
		
		A(Range::all(), 0, Range::all()) = A2d(Range::all(), Range::all());
		
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
	A = A/(DP(N[1]*N[2]*N[3]));  
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


