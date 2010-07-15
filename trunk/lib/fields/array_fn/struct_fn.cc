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
/*! \file struct_fn.cc
 * 
 * @sa struct_fn.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 * @bug  No known bugs
 */ 
 
#include "struct_fn.h"


//*********************************************************************************************



void Compute_dr
(
	string basis_type,
	int N[], 
	TinyVector<int,3> j1, TinyVector<int,3> j2, 
	TinyVector<DP,3> dr,
	DP xfactor[]
) 
{
	if (basis_type == "FOUR")					// Periodic dirn
	{
		if ((j2(0)-j1(0)) <= N[1]/2)
			dr(0) = (j2(0)-j1(0)) * xfactor[1];
			
		else
			dr(0) = (j2(0)-j1(0)-N[1]) * xfactor[1];
	}
	
	else if (basis_type == "SCFT")
		dr(0) = (j2(0)-j1(0)) * xfactor[1];		// Nonperiodic
		
	
	// Y, Z direction- periodic for both FOUR & SCFT
	if ((basis_type == "FOUR") || (basis_type == "SCFT"))
	{			
		if ((j2(1)-j1(1)) <= N[2]/2)
			dr(1) = (j2(1)-j1(1)) * xfactor[2];
			
		else
			dr(1) = (j2(1)-j1(1)-N[2]) * xfactor[2];
			
		//
		if ((j2(2)-j1(2)) <= N[3]/2)
			dr(2) = (j2(2)-j1(2)) * xfactor[3];
			
		else
			dr(2) = (j2(2)-j1(2)-N[3]) * xfactor[3];
	}
}



//*********************************************************************************************

void Compute_dSt
(
	string basis_type,
	int N[],
	TinyVector<int,3> j1, TinyVector<int,3> j2,
	TinyVector<DP,3> V1, TinyVector<DP,3> V2,
	int	q_min, int q_max,
	Array<DP,3> St, 
	Array<int,1> count,
	DP xfactor[]
)
{
	
	TinyVector<DP,3>	dr, dr_hat;
	TinyVector<DP,3>	dV;	
	DP					abs_dr, dV_pll, dV_perp;
	
	
	Compute_dr(basis_type, N, j1, j2, dr, xfactor);
	
	abs_dr = sqrt(dot(dr,dr));
	dr_hat = dr(0)/abs_dr, dr(1)/abs_dr, dr(2)/abs_dr;
	
	dV = V2(0)-V1(0), V2(1)-V1(1), V2(2)-V1(2);
	
	dV_pll = dot(dV, dr_hat);
	dV_perp = sqrt(dot(dV,dV) - pow2(dV_pll));
	

	int  rindex = (int) ceil(abs_dr);
	
	// Skip when P and Q overlap, 
	// or when the distance is larger than max_radius_inside_real_space.
	if ((abs_dr > MYEPS) && (rindex < (St.length)(0)))
		for (int q=q_min; q<=q_max; q++)
		{
			St(rindex,q,0) += pow(dV_pll, q);
			St(rindex,q,1) += pow(dV_perp, q);
			count(rindex) = count(rindex)+1;
		}	
	
}



//*********************************************************************************************

void Compute_dSt
(
	string basis_type,
	int N[],
	TinyVector<int,3> j1, TinyVector<int,3> j2,
	DP F1, DP F2,
	int	q_min, int q_max,
	Array<DP,2> St, 
	Array<int,1> count,
	DP xfactor[]
)
{	
	TinyVector<DP,3>	dr;	
	DP					abs_dr;
	
	Compute_dr(basis_type, N, j1, j2, dr, xfactor);
	abs_dr = sqrt(dot(dr,dr));
	
	int rindex = (int) ceil(abs_dr);
	
	// Skip when P and Q overlap, 
	// or when the distance is larger than max_radius_inside_real_space.
	if ((abs_dr > MYEPS) && (rindex < (St.length)(0)))
		for (int q=q_min; q<=q_max; q++)
		{
			St(rindex,q) += pow((F2-F1), q);
			count(rindex) = count(rindex)+1;
		}	
	
}


//******************************************************************************

void Compute_structure_function
(
	string basis_type,
	int N[], 
	Array<complx,3> Vx, Array<complx,3> Vy, Array<complx,3> Vz,
	int q_min, int q_max, 
	Array<DP,3> St, 
	DP xfactor[]
)
{	
}


//*********************************************************************************************

void Compute_structure_function
(
	string basis_type,
	int N[], 
	Array<complx,3> F,
	int q_min, int q_max, 
	Array<DP,2> St, 
	DP xfactor[]
)
{
	
}




//*********************************   End of struct_fn.cc **********************


