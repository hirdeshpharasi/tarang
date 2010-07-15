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
/*! \file planar_struct_fn.cc
 * 
 * @sa planar_struct_fn.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 * @bug  No known bugs
 */ 
 
#include "planar_struct_fn.h"
#include "struct_fn.h"


/**********************************************************************************************

	Planar structure function S(z, r, q) if anisotropy along z.
	Similar scheme for other anisotropy dirn.

***********************************************************************************************/

void Compute_dSt_planar
(
	string basis_type,
	int N[],
	TinyVector<int,3> j1, TinyVector<int,3> j2,
	TinyVector<DP,3> V1, TinyVector<DP,3> V2,
	int	q_min, int q_max,
	int kpll,
	Array<DP,4> St_planar, 
	Array<int,2> count,
	DP xfactor[]
)
{	
	
}

//*********************************************************************************************

void Compute_dSt_planar
(
	string basis_type,
	int N[],
	TinyVector<int,3> j1, TinyVector<int,3> j2,
	DP F1, DP F2,
	int	q_min, int q_max,
	int kpll,
	Array<DP,3> St_planar, 
	Array<int,2> count,
	DP xfactor[]
)
{	
}


//*********************************************************************************************


void Compute_planar_structure_function
(
	string basis_type,
	int N[], 
	Array<complx,3> Vx, Array<complx,3> Vy, Array<complx,3> Vz,
	int q_min, int q_max, 
	Array<DP,4> St_planar, 
	DP xfactor[]
)
{

}



//*********************************************************************************************


void Compute_planar_structure_function
(
	string basis_type,
	int N[], 
	Array<complx,3> F,
	int q_min, int q_max, 
	Array<DP,3> St_planar, 
	DP xfactor[]
)
{
}


//*********************************   End of struct_fn.cc *************************************




