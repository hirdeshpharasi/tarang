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
/*! \file  fill_ring.cc
 * 
 * @brief  Fill rings in a region.
 * 
 * We have two types of rings: (a) from a spherical shell, (b) from a cylindrical shell. <BR>
 *
 * (a) Spherical ring: ring_shell_index and sector_index are passed on as parameters.  
 *		The function picks inner-radius, outer-radius, left-angle and right-angle. <BR>
 *		Ring(n,m) = \f$ (R(n-1),\Theta(m-1); \R(n), \Theta(m)] \f$; (n,m>0) <BR>
 *		Exclude the modes in the inner sphere and include the modes in the outer sphere. <BR>
 *		Exclude the modes on the lef-angle surface, but include the modes on the right-angle
 *		surface.  Exception: \f$ Theta=0 \f$ modes however are kept in the first sector.
 *
 *	(b) Cylindrical ring: cylinder_shell_index and slab_index are passed on as parameters.
 *		The function picks inner-radius, outer-radius, h_lower, and h_upper. <BR>
 *		Cyl(n,m) = \f$ (R(n-1), H(m-1); R(n), H(m)] \f$; (n,m>0) <BR>
 *		Exclude the modes in the inner sphere and include the modes in the outer sphere. <BR>
 *		Exclude the modes on the lower surface (h_lower) and include the modes on the 
 *		upper surface.  Exception modes of H = lowest_h surface are kept in the first slab.
 *
 * @sa void Fill_array_ring(string basis_type, int N[], Array<complx,2> A, Array<complx,2> B,
 *			DP inner_radius, DP outer_radius, DP left_angle, DP right_angle, DP kfactor[])	
 *
 * @sa void Fill_array_cylinder_ring(string basis_type, int N[], 
 *			Array<complx,3> A, Array<complx,3> B, DP inner_radius, DP outer_radius, 
 *			DP h_lower, DP h_upper, DP kfactor[])	
 *
 * @sa universal_ET.cc
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug	 No known bugs
 */

#include "../IncVF.h"
#include "../IncSF.h"


//*********************************************************************************************

// Spherical rings
// Fluid
void IncVF::Fill_ring(int ring_shell_from_index, int sector_from_index)						
{	
	DP inner_radius = (*ring_shell_radius)(ring_shell_from_index - 1);
	DP outer_radius = (*ring_shell_radius)(ring_shell_from_index);

	
	DP left_angle = (*sector_angle_ring_tr)(sector_from_index-1);
	DP right_angle = (*sector_angle_ring_tr)(sector_from_index);

	
	Fill_array_ring(basis_type, N, *V1, *G1, inner_radius, outer_radius, 
						left_angle, right_angle, kfactor);
						
	Fill_array_ring(basis_type, N, *V2, *G2, inner_radius, outer_radius, 
						left_angle, right_angle, kfactor);
						
	Fill_array_ring(basis_type, N, *V3, *G3, inner_radius, outer_radius, 
						left_angle, right_angle, kfactor);
}


						
// Passive scalar

void IncVF::Fill_ring(int ring_shell_from_index, int sector_from_index, IncSF& T)						
{
	DP inner_radius = (*ring_shell_radius)(ring_shell_from_index - 1);
	DP outer_radius = (*ring_shell_radius)(ring_shell_from_index);

	DP left_angle = (*sector_angle_ring_tr)(sector_from_index-1);
	DP right_angle = (*sector_angle_ring_tr)(sector_from_index);
				
	Fill_array_ring(basis_type, N, *T.F, *G1,  inner_radius, outer_radius, 
						left_angle, right_angle, kfactor);	
}



// MHD B field

void IncVF::Fill_ring(int ring_shell_from_index, int sector_from_index, IncVF& W)						
{
	DP inner_radius = (*ring_shell_radius)(ring_shell_from_index - 1);
	DP outer_radius = (*ring_shell_radius)(ring_shell_from_index);
		
	DP left_angle = (*sector_angle_ring_tr)(sector_from_index-1);
	DP right_angle = (*sector_angle_ring_tr)(sector_from_index);	
	

	Fill_array_ring(basis_type, N, *W.V1, *G1, inner_radius, outer_radius, 
						left_angle, right_angle, kfactor);
						
	Fill_array_ring(basis_type, N, *W.V2, *G2, inner_radius, outer_radius, 
						left_angle, right_angle, kfactor);
						
	Fill_array_ring(basis_type, N, *W.V3, *G3, inner_radius, outer_radius, 
						left_angle, right_angle, kfactor);
}

//*********************************************************************************************

// Cylindrical rings
// Fluid

void IncVF::Fill_cylinder_ring(int cylinder_shell_from_index, int slab_from_index)						
{	

	DP inner_radius = (*cylinder_shell_radius)(cylinder_shell_from_index - 1);
	DP outer_radius = (*cylinder_shell_radius)(cylinder_shell_from_index);

	
	DP h_lower = (*cylinder_kpll_array_tr)(slab_from_index-1);
	DP h_upper = (*cylinder_kpll_array_tr)(slab_from_index);
	
	
	Fill_array_cylinder_ring(basis_type, N, *V1, *G1, inner_radius, outer_radius, 
						h_lower, h_upper, kfactor);
						
	Fill_array_cylinder_ring(basis_type, N, *V2, *G2, inner_radius, outer_radius, 
						h_lower, h_upper, kfactor);
						
	Fill_array_cylinder_ring(basis_type, N, *V3, *G3, inner_radius, outer_radius, 
						h_lower, h_upper, kfactor);
}


// Scalar

void IncVF::Fill_cylinder_ring(int cylinder_shell_from_index, int slab_from_index, IncSF& T)						
{	

	DP inner_radius = (*cylinder_shell_radius)(cylinder_shell_from_index - 1);
	DP outer_radius = (*cylinder_shell_radius)(cylinder_shell_from_index);

	
	DP h_lower = (*cylinder_kpll_array_tr)(slab_from_index-1);
	DP h_upper = (*cylinder_kpll_array_tr)(slab_from_index);
	
	
	Fill_array_cylinder_ring(basis_type, N, *T.F, *G1,  inner_radius, outer_radius, 
							h_lower, h_upper, kfactor);
}

// MHD B field

void IncVF::Fill_cylinder_ring(int cylinder_shell_from_index, int slab_from_index, IncVF& W)						
{	

	DP inner_radius = (*cylinder_shell_radius)(cylinder_shell_from_index - 1);
	DP outer_radius = (*cylinder_shell_radius)(cylinder_shell_from_index);

	
	DP h_lower = (*cylinder_kpll_array_tr)(slab_from_index-1);
	DP h_upper = (*cylinder_kpll_array_tr)(slab_from_index);
	
	
	Fill_array_cylinder_ring(basis_type, N, *W.V1, *G1, inner_radius, outer_radius, 
						h_lower, h_upper, kfactor);
						
	Fill_array_cylinder_ring(basis_type, N, *W.V2, *G2, inner_radius, outer_radius, 
						h_lower, h_upper, kfactor);
						
	Fill_array_cylinder_ring(basis_type, N, *W.V3, *G3, inner_radius, outer_radius, 
						h_lower, h_upper, kfactor);
}

//********************************  End of fill_rings.cc **************************************





