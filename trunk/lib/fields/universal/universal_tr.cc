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


/*! \file  universal_tr.cc
 * 
 * @sa	universal_fn.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 @ @bug  No known bug
 */
 
#include "universal_tr.h"


//*********************************************************************************************
//
//     Init fftw plans
//
//*********************************************************************************************


void Init_fftw_plan_array(string basis_type, int N[], Array<complx,3> A)
{
	
	if (basis_type == "FOUR")
		Init_fftw_plan_FOUR(N, A);
	
	else if (basis_type == "SCFT")
		Init_fftw_plan_SCFT(N, A);
}


//*********************************************************************************************
//
// Forward transform
//
//*********************************************************************************************


void Forward_transform_array
(
	string basis_type, 
	int N[],  
	Array<complx,3> A, 
	Array<complx,3> temp_r, 
	int parity
)										
{	
	if (basis_type == "FOUR")
		ArrayFFT_FOUR(N, A, temp_r);
		
	else if (basis_type == "SCFT")
	{
		if (parity == 0)					// even parity
			ArrayCFT_SCFT(costr_plan_SCFT, r2c_plan_SCFT, r2c_1d_plan_SCFT, N, A, temp_r);
			
		else if (parity == 1)				// odd parity
			ArraySFT_SCFT(sintr_plan_SCFT, r2c_plan_SCFT, r2c_1d_plan_SCFT, N, A, temp_r);
	}		
}


//
// Forward transform_transpose_order A: T(Atr) -> A
//			


void Forward_transform_array_transpose_order
(
	string basis_type, 
	int N[],
	Array<complx,3> Atr, 
	Array<complx,3> A, 
	int parity
)
{
	if (basis_type == "FOUR")
		ArrayFFT_FOUR_transpose_order(N, Atr, A);
		
	else if (basis_type == "SCFT") 
	{
		if (parity == 0)					// even parity
			ArrayCFT_SCFT_transpose_order(costr_plan_SCFT, r2c_plan_SCFT, r2c_1d_plan_SCFT, 
													N, Atr, A);
			
		else if (parity == 1)				// odd parity
			ArraySFT_SCFT_transpose_order(sintr_plan_SCFT, r2c_plan_SCFT, r2c_1d_plan_SCFT, 
													N, Atr, A);
	}		
}



//*********************************************************************************************
//
// Inverse transform
//
//*********************************************************************************************

void Inverse_transform_array
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> temp_r, 
	int parity
) 							
{
	if (basis_type == "FOUR")
		ArrayIFFT_FOUR(N, A, temp_r);
	
	else if (basis_type == "SCFT") 
	{
		if (parity == 0)					// even parity
			ArrayICFT_SCFT(icostr_plan_SCFT, c2r_plan_SCFT, c2r_1d_plan_SCFT, N, A, temp_r);
			
		else if (parity == 1)				// odd parity
			ArrayISFT_SCFT(isintr_plan_SCFT, c2r_plan_SCFT, c2r_1d_plan_SCFT, N, A, temp_r);
	}
}


//
//  Inverse transform_transpose_order: I(A) -> Atr
//			

void Inverse_transform_array_transpose_order
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> Atr, 
	int parity
)
{
	if (basis_type == "FOUR")
		ArrayIFFT_FOUR_transpose_order(N, A, Atr);
	
	else if (basis_type == "SCFT") 
	{
		if (parity == 0)					// even parity
			ArrayICFT_SCFT_transpose_order(icostr_plan_SCFT, c2r_plan_SCFT, c2r_1d_plan_SCFT, 
												N, A, Atr);
			
		else if (parity == 1)				// odd parity
			ArrayISFT_SCFT_transpose_order(isintr_plan_SCFT, c2r_plan_SCFT, c2r_1d_plan_SCFT, 
												N, A, Atr);
	}		
}



//*********************************************************************************************
//
//     		Derivative[A] -> B
//	
//*********************************************************************************************


void Xderiv
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> B, 
	DP kfactor[], 
	int parity
)
{
	if (basis_type == "FOUR")
		Xderiv_FOUR(N, A, B, kfactor);
		
	else if (basis_type == "SCFT") 
	{
		if (parity == 0)					// even parity
			Xderiv_Cos_SCFT(N, A, B, kfactor);
			
		else if (parity == 1)				// odd parity
			Xderiv_Sin_SCFT(N, A, B, kfactor);
	}		
}	

// Yderive:    B(k) =	i Ky x A(k)


void Yderiv
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> B, 
	DP kfactor[], 
	int parity
)
{
	if (basis_type == "FOUR")
		Yderiv_FOUR(N, A, B, kfactor);
		
	else if (basis_type == "SCFT")
		Yderiv_SCFT(N, A, B, kfactor);
}


// Zderive	 B(k) =	iKz x A(k)


void Zderiv
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	Array<complx,3> B, 
	DP kfactor[], 
	int parity
)
{
	if (basis_type == "FOUR")
		Zderiv_FOUR(N, A, B, kfactor);
		
	else if (basis_type == "SCFT")
		Zderiv_SCFT(N, A, B, kfactor);
}

//********************************  End of universal_fn.cc   **********************************



