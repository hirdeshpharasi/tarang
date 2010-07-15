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


/*! \file  universal_ET.cc
 * 
 * @sa	universal_ET.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bug	  No known bugs
 */
 
#include "universal_ET.h"	 


//*********************************************************************************************


void Fill_array_shell
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP kfactor[]
)					
{	
	DP kmag;
	B = 0.0;
	for (int l1=0;  l1<local_N1; l1++)
		for (int l2=0;  l2<N[2]; l2++) 
			for (int l3=0;  l3<=N[3]/2; l3++) 
			{
				kmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
				
				if ( (kmag > inner_radius) && (kmag <= outer_radius) ) 
					(B)(l1,l2,l3) = (A)(l1,l2,l3); 
			}
}

//
//

void Fill_array_ring
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP left_angle, DP right_angle, 
	DP kfactor[]
)						
{	
	DP kkmag, theta;
	
	B = 0.0;
	for (int l1=0;  l1<local_N1; l1++)
		for (int l2=0;  l2<N[2]; l2++) 
			for (int l3=0;  l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
				
				theta = AnisKvect_polar_angle(basis_type, l1, l2, l3, N, kfactor);
				
				if ( (kkmag > inner_radius) && (kkmag <= outer_radius) )
					if ( ((theta > left_angle) && (theta <= right_angle)) 
								|| ((abs(theta) < MYEPS) && (left_angle < MYEPS)) )
						(B)(l1,l2,l3) = (A)(l1,l2,l3); 
			}
}

//
//

void Fill_array_cylinder_ring
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP h_lower, DP h_upper, 
	DP kfactor[]
)						
{	
	DP kkpll, kkperp;
	
	B = 0.0;
	for (int l1=0;  l1<N[1]; l1++)
		for (int l2=0;  l2<N[2]; l2++) 
			for (int l3=0;  l3<=N[3]/2; l3++) 
			{
				kkpll = AnisKpll(basis_type, l1, l2, l3, N, kfactor);
				
				kkperp = AnisKperp(basis_type, l1, l2, l3, N, kfactor);
				
				if ( (kkperp > inner_radius) && (kkperp <= outer_radius) )
					if ( ((kkpll > h_lower) && (kkpll <= h_upper)) 
								|| ((abs(kkpll) < MYEPS) && (h_lower < MYEPS)) )
						(B)(l1,l2,l3) = (A)(l1,l2,l3); 
			}
}


//*********************************************************************************************
	


DP	Shell_mult_single
(
	string basis_type,
	string alias_switch, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		return Shell_mult_single_FOUR(alias_switch, N, A, B, 
										inner_radius, outer_radius, kfactor);
		
	else if (basis_type == "SCFT")
		return Shell_mult_single_SCFT(alias_switch, N, A, B, 
										inner_radius, outer_radius, kfactor);
		
	else
		return 0;  // for -Wall		
}

//
//
	

void Shell_mult_single_real_imag
(
	string basis_type, 
	string alias_switch,
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP& total_real, DP& total_imag,
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		Shell_mult_single_real_imag_FOUR(alias_switch, N, A, B, inner_radius, outer_radius, 
													total_real, total_imag, kfactor);
	else if (basis_type == "SCFT")
		Shell_mult_single_real_imag_SCFT(alias_switch, N, A, B, inner_radius, outer_radius, 
													total_real, total_imag, kfactor);
}		


//*********************************************************************************************

void Shell_mult_all
(
	string basis_type, 
	string alias_switch,
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		Shell_mult_all_FOUR(alias_switch, N, A, B, shell_radius_array, result, kfactor);
		
	else if (basis_type == "SCFT")
		Shell_mult_all_SCFT(alias_switch, N, A, B, shell_radius_array, result, kfactor);
}


void Shell_mult_all_real_imag
(
	string basis_type, 
	string alias_switch,
	int N[], 
	Array<complx,3> A,  Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result_real, Array<DP,1> result_imag, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		Shell_mult_all_real_imag_FOUR(alias_switch, N, A, B, shell_radius_array, 
										result_real, result_imag, kfactor);
										
	else if (basis_type == "SCFT")
		Shell_mult_all_real_imag_SCFT(alias_switch, N, A, B, shell_radius_array, 
										result_real, result_imag, kfactor);
}		

// 
// DOT PROD
//

void Shell_mult_all
(
	string basis_type, 
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		Shell_mult_all_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										shell_radius_array, result, kfactor);
		
	else if (basis_type == "SCFT")
		Shell_mult_all_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										shell_radius_array, result, kfactor);
}


void Shell_mult_all_real_imag
(
	string basis_type, 
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result_real, Array<DP,1> result_imag, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		Shell_mult_all_real_imag_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										shell_radius_array, result_real, result_imag, kfactor);
										
	else if (basis_type == "SCFT")
		Shell_mult_all_real_imag_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										shell_radius_array, result_real, result_imag, kfactor);
}		


//********************************************************************************************* 


void Ring_mult_all
(
	string basis_type, 
	string alias_switch,
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, Array<DP, 1> sector_angle_array, 
	Array<DP,2> result, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		Ring_mult_all_FOUR(alias_switch, N, A, B, shell_radius_array, sector_angle_array, 
										result, kfactor);
		
	else if (basis_type == "SCFT")
		Ring_mult_all_SCFT(alias_switch, N, A, B, shell_radius_array, sector_angle_array, 
										result, kfactor);
}


void Ring_mult_all_real_imag
(
	string basis_type, 
	string alias_switch,
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, Array<DP, 1> sector_angle_array,
	Array<DP,2> result_real, Array<DP,2> result_imag, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		Ring_mult_all_real_imag_FOUR(alias_switch, N, A, B, shell_radius_array, 
									sector_angle_array,result_real, result_imag, kfactor);
									
	else if (basis_type == "SCFT")
		Ring_mult_all_real_imag_SCFT(alias_switch, N, A, B, shell_radius_array, 
									sector_angle_array, result_real, result_imag, kfactor);
}	

//
// DOT PROD
//

void Ring_mult_all
(
	string basis_type, 
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, Array<DP, 1> sector_angle_array, 
	Array<DP,2> result, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		Ring_mult_all_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, shell_radius_array, 
							sector_angle_array, result, kfactor);
		
	else if (basis_type == "SCFT")
		Ring_mult_all_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, shell_radius_array, 
							sector_angle_array, result, kfactor);
}

//
//

void Ring_mult_all_real_imag
(
	string basis_type, 
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, Array<DP, 1> sector_angle_array,
	Array<DP,2> result_real, Array<DP,2> result_imag, 
	DP kfactor[]
)
{
	if (basis_type == "FOUR") 
		Ring_mult_all_real_imag_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										shell_radius_array,  sector_angle_array, 
										result_real, result_imag, kfactor);
								
	else if (basis_type == "SCFT")
		Ring_mult_all_real_imag_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										shell_radius_array,  sector_angle_array, 
										result_real, result_imag, kfactor);
}		

//*********************************************************************************************

void Cyl_ring_mult_all
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> cylinder_shell_radius_array,  Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> result, 
	DP kfactor[]
)
{

	if (basis_type == "FOUR") 
		Cyl_ring_mult_all_FOUR(alias_switch, N, A, B, cylinder_shell_radius_array,  
								cylinder_kpll_array, result, kfactor);
		
	else if (basis_type == "SCFT")
		Cyl_ring_mult_all_SCFT(alias_switch, N, A, B, cylinder_shell_radius_array, 
								cylinder_kpll_array, result, kfactor);
}		

//*********************************************************************************************

void Cyl_ring_mult_all_real_imag
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> cylinder_shell_radius_array, Array<DP, 1> cylinder_kpll_array,
	Array<DP,2> result_real, Array<DP,2> result_imag, 
	DP kfactor[]
)
{

	if (basis_type == "FOUR") 
		Cyl_ring_mult_all_real_imag_FOUR(alias_switch, N, A, B, cylinder_shell_radius_array,  
								cylinder_kpll_array, result_real, result_imag, kfactor);
		
	else if (basis_type == "SCFT")
		Cyl_ring_mult_all_real_imag_SCFT(alias_switch, N, A, B, cylinder_shell_radius_array, 
								cylinder_kpll_array, result_real, result_imag, kfactor);
}	

//*********************************************************************************************

void Cyl_ring_mult_all
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array,  Array<DP, 1> cylinder_kpll_array,
	Array<DP,2> result, 
	DP kfactor[]
)
{

	if (basis_type == "FOUR") 
		Cyl_ring_mult_all_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz,		
						cylinder_shell_radius_array, cylinder_kpll_array, result, kfactor);
		
	else if (basis_type == "SCFT")
		Cyl_ring_mult_all_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
						cylinder_shell_radius_array, cylinder_kpll_array, result, kfactor);
}	


//*********************************************************************************************

void Cyl_ring_mult_all_real_imag
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array, Array<DP, 1> cylinder_kpll_array,
	Array<DP,2> result_real, Array<DP,2> result_imag, 
	DP kfactor[]
)
{

	if (basis_type == "FOUR") 
		Cyl_ring_mult_all_real_imag_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
								cylinder_shell_radius_array, cylinder_kpll_array, 
								result_real, result_imag, kfactor);
		
	else if (basis_type == "SCFT")
		Cyl_ring_mult_all_real_imag_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
								cylinder_shell_radius_array, cylinder_kpll_array, 
								result_real, result_imag, kfactor);
}		


//*********************************************************************************************

void Shell_mult_all_imagVW_B0
(
	string basis_type,
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result, 
	DP kfactor[]
)
{

	if (basis_type == "FOUR") 
		Shell_mult_all_imagVW_B0_FOUR(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz,
										shell_radius_array, result, kfactor); 
		
	else if (basis_type == "SCFT")
		Shell_mult_all_imagVW_B0_SCFT(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz,
										shell_radius_array, result, kfactor);
}

//
//

void Ring_mult_all_imagVW_B0
(
	string basis_type,
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> ring_shell_radius_array,  Array<DP, 1> sector_angle_array, 
	Array<DP,2> result, 
	DP kfactor[]
)
{

	if (basis_type == "FOUR") 
		Ring_mult_all_imagVW_B0_FOUR(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz,
										ring_shell_radius_array, sector_angle_array,
										result, kfactor); 
		
	else if (basis_type == "SCFT")
		Ring_mult_all_imagVW_B0_SCFT(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz,
										ring_shell_radius_array, sector_angle_array,
										result, kfactor);
}

//
//

void Cyl_ring_mult_all_imagVW_B0
(
	string basis_type,
	string alias_switch, 
	int N[], 
	TinyVector<DP,3> B0,
	Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
	Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
	Array<DP, 1> cylinder_shell_radius_array,  Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> result, 
	DP kfactor[]
)
{

	if (basis_type == "FOUR") 
		Cyl_ring_mult_all_imagVW_B0_FOUR(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz,
										cylinder_shell_radius_array, cylinder_kpll_array,
										result, kfactor); 
		
	else if (basis_type == "SCFT")
		Cyl_ring_mult_all_imagVW_B0_SCFT(alias_switch, N, B0, Ax, Ay, Az, Bx, By, Bz,
										cylinder_shell_radius_array, cylinder_kpll_array,
										result, kfactor); 
}

//*********************************************************************************************

DP Shell_mult_vorticity
(
 string basis_type,
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 DP inner_radius, DP outer_radius, 
 DP kfactor[]
 )
{
	if (basis_type == "FOUR") 
		return Shell_mult_vorticity_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										 inner_radius, outer_radius, kfactor);
	
	else if (basis_type == "SCFT")
		return Shell_mult_vorticity_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
										 inner_radius, outer_radius, kfactor);
	
	else
		return 0;  // for -Wall	
	
}


DP Shell_mult_vector_potential
(
 string basis_type, 
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
 DP inner_radius, DP outer_radius, 
 DP kfactor[]
 )
{
	if (basis_type == "FOUR") 
		return Shell_mult_vector_potential_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
												inner_radius, outer_radius, kfactor);
	
	else if (basis_type == "SCFT")
		return Shell_mult_vector_potential_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
												inner_radius, outer_radius, kfactor);
	
	else
		return 0;  // for -Wall	
	
}




void Shell_mult_vorticity_all
(
 string basis_type,
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 Array<DP, 1> shell_radius_array, Array<DP,1> result, 
 DP kfactor[]
 )
{
	if (basis_type == "FOUR") 
		Shell_mult_vorticity_all_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
									  shell_radius_array, result, kfactor);
	
	else if (basis_type == "SCFT")
		Shell_mult_vorticity_all_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
									  shell_radius_array, result, kfactor);
	
}


void Shell_mult_vector_potential_all
(
 string basis_type, 
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
 Array<DP, 1> shell_radius_array, Array<DP,1> result, 
 DP kfactor[]
 )
{
	if (basis_type == "FOUR") 
		Shell_mult_vector_potential_all_FOUR(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
											 shell_radius_array, result, kfactor);
	
	else if (basis_type == "SCFT")
		Shell_mult_vector_potential_all_SCFT(alias_switch, N, Ax, Ay, Az, Bx, By, Bz, 
												shell_radius_array, result, kfactor);
	
}




//***************************   End of universal_ET.cc   **************************************


