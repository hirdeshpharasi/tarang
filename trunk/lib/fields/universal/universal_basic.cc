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


/*! \file  universal_basic.cc
 * 
 * @sa	universal_basic.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 * @bug	 No known bugs
 */
 
 
#include "universal_basic.h"


//*********************************************************************************************
 
// No of modes in the shell[inner_radius, outer_radius).
int Get_number_modes_in_shell
(
	string basis_type, 
	int N[], 
	DP inner_radius, DP outer_radius, 
	DP kfactor[]
)
{

	if (basis_type == "FOUR") 
		return  Get_number_modes_in_shell_FOUR(N, inner_radius, outer_radius, kfactor);
		
	else if (basis_type == "SCFT")
		return Get_number_modes_in_shell_SCFT(N, inner_radius, outer_radius, kfactor);
		
	else
		return 0;		// for -Wall

}

/// Grid wavenumber at (l1, l2, l3).
void Wavenumber
(
	string basis_type, 
	int l1, int l2, int l3, 
	int N[], 
	DP kfactor[], 
	TinyVector<DP,3> &kk
)
{
	if (basis_type == "FOUR") 
		Wavenumber_FOUR(l1, l2, l3, N, kfactor, kk);
		
	else if (basis_type == "SCFT")
		Wavenumber_SCFT(l1, l2, l3, N, kfactor, kk);
}


// Get Modal helicity for (l1,l2,l3).
DP Get_Modal_helicity
(
	string basis_type, 
	int l1, int l2, int l3, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		return Get_Modal_helicity_FOUR(l1, l2, l3, N, Ax, Ay, Az, kfactor);
		
	else if (basis_type == "SCFT")
		return Get_Modal_helicity_SCFT(l1, l2, l3, N, Ax, Ay, Az, kfactor);
		
	else
		return 0;	// Safe exit 	
}									


// Computes vorticity of mode with grid index (l1, l2, l3).
void Compute_Modal_vorticity
(
	string basis_type, 
	int l1, int l2, int l3, 
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
	DP kfactor[], 
	TinyVector<complx,3> &vorticity
)
{
	if (basis_type == "FOUR") 
		Compute_Modal_vorticity_FOUR(l1, l2, l3, N, Ax, Ay, Az, kfactor, vorticity);
		
	else if (basis_type == "SCFT")
		Compute_Modal_vorticity_SCFT(l1, l2, l3, N, Ax, Ay, Az, kfactor, vorticity);
}	

// Computes vorticity of mode with grid index (l1, l2, l3).
void Compute_Modal_vorticity_y_component
(
 string basis_type, 
 int l1, int l2, int l3, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az, 
 DP kfactor[], 
 complx &vort_y
 )
{
	if (basis_type == "FOUR") 
		Compute_Modal_vorticity_y_component_FOUR(l1, l2, l3, N, Ax, Ay, Az, kfactor, vort_y);
	
	else if (basis_type == "SCFT")
		Compute_Modal_vorticity_y_component_SCFT(l1, l2, l3, N, Ax, Ay, Az, kfactor, vort_y);
}	


//*********************************************************************************************


void Array_mult_ksqr(string basis_type, int N[], Array<complx,3> A, DP kfactor[])
{
	if (basis_type == "FOUR")
		Array_mult_ksqr_FOUR(N, A, kfactor);
		
	else if (basis_type == "SCFT") 
		Array_mult_ksqr_SCFT(N, A, kfactor);
		
}


void Array_divide_ksqr(string basis_type, int N[], Array<complx,3> A, DP kfactor[])
{
	if (basis_type == "FOUR")
		Array_divide_ksqr_FOUR(N, A, kfactor);
		
	else if (basis_type == "SCFT") 
		Array_divide_ksqr_SCFT(N, A, kfactor);
}



void Array_exp_ksqr(string basis_type, int N[], Array<complx,3> A, DP factor, DP kfactor[])
{
	if (basis_type == "FOUR")
		Array_exp_ksqr_FOUR(N, A, factor, kfactor);
		
	else if (basis_type == "SCFT") 
		Array_exp_ksqr_SCFT(N, A, factor, kfactor);
}


void Array_exp_ksqr
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	DP factor, DP hyper_factor,
	 DP kfactor[]
)
{
	if (basis_type == "FOUR")
		Array_exp_ksqr_FOUR(N, A, factor, hyper_factor, kfactor);
		
	else if (basis_type == "SCFT") 
		Array_exp_ksqr_SCFT(N, A, factor, hyper_factor, kfactor);
}

//*********************************************************************************************

void Array_mult_V0_khat_sqr
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, 
	TinyVector<DP,3> V0, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR")
		Array_mult_V0_khat_sqr_FOUR(N, A, V0, kfactor);
		
	else if (basis_type == "SCFT") 
		Array_mult_V0_khat_sqr_SCFT(N, A, V0, kfactor);
}


//********************************  End of universal_basis.cc   *******************************



