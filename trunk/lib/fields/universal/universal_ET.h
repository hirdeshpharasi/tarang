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


/*! \file  universal_ET.h
 * 
 * @brief  Universal functions based on basis fns for Energy transfers
 *
 *  Shell & ring functions for energy transfer calculations.
 * 
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bug  No known bugs
 */


#ifndef _UNIVERAL_ET_H
#define _UNIVERAL_ET_H

#include <mpi.h>
#include <fftw3-mpi.h>

#include <blitz/array.h>
#include <complex>
#include <cmath>
#include <string>

#ifdef BZ_HAVE_STD
#include <fstream>
#else
#include <fstream.h>
#endif


#include "../../four/four_ET.h"
#include "../../scft/scft_ET.h"

#include "universal_inline.h"


using namespace blitz;

//*********************************************************************************************	

	
///  3D: B(k) = A(k)  for k = [inner_radius, outer_radius); rest 0.
void Fill_array_shell
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP kfactor[]
);

///  3D: B(k) = A(k)  for k = [inner_radius, inner_angle; outer_radius, outer_angle); rest 0.
/// The fist sector also contains \f$ \theta = 0 \f$.
void Fill_array_ring
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP left_angle, DP right_angle, 
	DP kfactor[]
);	
								


///  3D: B(k) = A(k)  for k = (inner_radius, h_lower; outer_radius, h_upper]; rest 0.
/// The fist slab also contains \f$ k_{||} = 0 \f$.	
void Fill_array_cylinder_ring
(
	string basis_type, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP h_lower, DP h_upper, 
	DP kfactor[]
);


//*********************************************************************************************

/// Compute product  Real(A.B*) with appropriate factor for energy for a 3D 
///		shell = (inner_radius, outer_radius].	
DP	Shell_mult_single
(
	string basis_type,
	string alias_switch, 
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	DP inner_radius, DP outer_radius, 
	DP kfactor[]
);
					

/// Compute product  Real(A.B*) with appropriate factor for energy for all 3D shells. 	
///		The result is in array 'result'.								
void Shell_mult_all
(
	string basis_type, 
	string alias_switch,
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, 
	Array<DP,1> result, 
	DP kfactor[]
);


														
								
/// Compute product  \f$ \Re(\vec{A} \cdot \vec{B}^{*}) \f$ with appropriate factor for energy 
///		for all 3D shells. 	
///		The result is in array 'result'.								
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
);

//*********************************************************************************************
								
/// Ring mult

/// Compute product  Real(A.B*) with appropriate factor for energy for all 3D rings. 	
///		The result is in array 'result'.																						
void Ring_mult_all
(
	string basis_type, 
	string alias_switch,
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> shell_radius_array, Array<DP, 1> sector_angle_array, 
	Array<DP,2> result, 
	DP kfactor[]
);
	

/// Compute product  \f$ \Re(\vec{A} \cdot \vec{B}^{*}) \f$ with appropriate factor for energy 
///		for all 3D rings. 	
///		The result is in array 'result'.								
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
);
								
							
//*********************************************************************************************

//
// For cylinder
//

void Cyl_ring_mult_all
(
	string basis_type,
	string alias_switch,
	int N[], 
	Array<complx,3> A, Array<complx,3> B, 
	Array<DP, 1> cylinder_shell_radius_array,  Array<DP, 1> cylinder_kpll_array, 
	Array<DP,2> result, 
	DP kfactor[]
);


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
);




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
);

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
);


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
);

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
 );


DP Shell_mult_vector_potential
(
 string basis_type, 
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
 DP inner_radius, DP outer_radius, 
 DP kfactor[]
 );


void Shell_mult_vorticity_all
(
 string basis_type,
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz,  
 Array<DP, 1> shell_radius_array, Array<DP,1> result, 
 DP kfactor[]
 );

void Shell_mult_vector_potential_all
(
 string basis_type, 
 string alias_switch, 
 int N[], 
 Array<complx,3> Ax, Array<complx,3> Ay, Array<complx,3> Az,
 Array<complx,3> Bx, Array<complx,3> By, Array<complx,3> Bz, 
 Array<DP, 1> shell_radius_array, Array<DP,1> result, 
 DP kfactor[]
 );

#endif


//***************************   End of universal_ET.cc   **************************************


