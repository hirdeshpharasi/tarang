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

/*! \file  fill_sphere_shell.cc
 * 
 * @brief  Fill sphere or shell in a region.
 * 
 * Sphere_index / shell_index is passed to the function as a parameter.  The function
 *		picks the inner and outer radii and fill the given region. <BR>
 *		Notation shell(n) = (R(n-1),R(n)] (excluding the modes in the inner sphere and 
 *		including the modes in the outer sphere. <BR>
 *		The last shell contains all the modes beyond the maximum inner sphere.
 *
 * @sa void Fill_array_shell(string basis_type, int N[], Array<complx,2> A, Array<complx,2> B,
 *		DP inner_radius, DP outer_radius, 	DP kfactor[])
 *
 * @sa universal_ET.cc
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "../IncVF.h"
#include "../IncSF.h"


//*********************************************************************************************

///  FOR FLUID --- INSIDE SPHERE FILLED
///	 Includes the modes of the outer surface.
void IncVF::Fill_in_sphere(int sphere_index)						
{
	DP radius = (*sphere_radius)(sphere_index);				// radius of the sphere

	Fill_array_shell(basis_type, N, *V1, *G1, 0, radius, kfactor);
	Fill_array_shell(basis_type, N, *V2, *G2, 0, radius, kfactor);
	Fill_array_shell(basis_type, N, *V3, *G3, 0, radius, kfactor);
}


// For Passive Scalar -- INSIDE SPHERE FILLED

void IncVF::Fill_in_sphere(int sphere_index, IncSF& T)						
{
		
	DP radius = (*sphere_radius)(sphere_index);				// radius of the sphere
			
	Fill_array_shell(basis_type, N, *T.F, *G1, 0, radius, kfactor);			
}


// For MHD -- INSIDE SPHERE FILLED

void IncVF::Fill_in_sphere(int sphere_index, IncVF& W)						
{
	DP radius = (*sphere_radius)(sphere_index);				// radius of the sphere
	
	Fill_array_shell(basis_type, N, *W.V1, *G1, 0, radius, kfactor);
	Fill_array_shell(basis_type, N, *W.V2, *G2, 0, radius, kfactor);
	Fill_array_shell(basis_type, N, *W.V3, *G3, 0, radius, kfactor);
}


/*====================================================================================

void IncVF::Fill_out_sphere(): -- OUTSIDE SPHERE FILLED
	 
// Modes on the surface are included outside the sphere

======================================================================================*/ 

// FOR FLUID --- OUTSIDE SPHERE FILLED

void IncVF::Fill_out_sphere(int sphere_index)						
{
	DP radius = (*sphere_radius)(sphere_index);				// radius of the sphere

	Fill_array_shell(basis_type, N, *V1, *G1, radius, INF_RADIUS, kfactor);
	Fill_array_shell(basis_type, N, *V2, *G2, radius, INF_RADIUS, kfactor);
	Fill_array_shell(basis_type, N, *V3, *G3, radius, INF_RADIUS, kfactor);
}


// For Passive Scalar 
void IncVF::Fill_out_sphere(int sphere_index, IncSF& T)						
{
	DP radius = (*sphere_radius)(sphere_index);				// radius of the sphere
			
	Fill_array_shell(basis_type, N, *T.F, *G1,  radius, INF_RADIUS, kfactor);
}


// For MHD 
void IncVF::Fill_out_sphere(int sphere_index, IncVF& W)						
{
	DP radius = (*sphere_radius)(sphere_index);				// radius of the sphere

	Fill_array_shell(basis_type, N, *W.V1, *G1, radius, INF_RADIUS, kfactor);
	Fill_array_shell(basis_type, N, *W.V2, *G2, radius, INF_RADIUS, kfactor);
	Fill_array_shell(basis_type, N, *W.V3, *G3, radius, INF_RADIUS, kfactor);
}



//*********************************************************************************************

/// Fluid -- fill shell
/// Include the modes of the outer surface, and exclude that of inner surface.
void IncVF::Fill_shell(int shell_from_index)						
{
	DP inner_radius = (*shell_radius)(shell_from_index - 1);
	DP outer_radius = (*shell_radius)(shell_from_index);
	
	Fill_array_shell(basis_type, N, *V1, *G1, inner_radius, outer_radius, kfactor);
	Fill_array_shell(basis_type, N, *V2, *G2, inner_radius, outer_radius, kfactor);
	Fill_array_shell(basis_type, N, *V3, *G3, inner_radius, outer_radius, kfactor);
}


// Passive scalar
void IncVF::Fill_shell(int shell_from_index, IncSF& T)						
{
	DP inner_radius = (*shell_radius)(shell_from_index - 1);
	DP outer_radius = (*shell_radius)(shell_from_index);
			
	Fill_array_shell(basis_type, N, *T.F, *G1,  inner_radius, outer_radius, kfactor);
}


// MHD B field
void IncVF::Fill_shell(int shell_from_index, IncVF& W)						
{
	DP inner_radius = (*shell_radius)(shell_from_index - 1);
	DP outer_radius = (*shell_radius)(shell_from_index);
	
	Fill_array_shell(basis_type, N, *W.V1, *G1, inner_radius, outer_radius, kfactor);
	Fill_array_shell(basis_type, N, *W.V2, *G2, inner_radius, outer_radius, kfactor);
	Fill_array_shell(basis_type, N, *W.V3, *G3, inner_radius, outer_radius, kfactor);
}


	//*********************************************************************************************

	///  FOR FLUID --- INSIDE SPHERE FILLED
	///	 For helicity flux calculations
void IncVF::Fill_in_sphere_vorticity(int sphere_index)						
{
	
	TinyVector<complx,3> vorticity;
	complx vort_y;
	
	DP radius = (*sphere_radius)(sphere_index);				// radius of the sphere
	
	DP kkmag;
	
	*G1 = 0.0;
	*G2 = 0.0;
	*G3 = 0.0;
	
	for (int l1=0;  l1<local_N1; l1++)
		for (int l2=0;  l2<N[2]; l2++) 
			for (int l3=0;  l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
				
				if ((kkmag > 0) && (kkmag <= radius))
				{
					if (N[2] > 1)
					{	
						Compute_Modal_vorticity(basis_type, l1, l2, l3, N, 
												*V1, *V2, *V3, kfactor, vorticity);
						
						(*G1)(l1, l2, l3) = vorticity(0);
						(*G2)(l1, l2, l3) = vorticity(1);
						(*G3)(l1, l2, l3) = vorticity(2);
					}	
					
					else if (N[2] ==1)
					{
						Compute_Modal_vorticity_y_component(basis_type, l1, 0, l3, N, 
															*V1, *V1, *V3, kfactor, vort_y);
						
						(*G1)(l1, 0, l3) = vort_y;						
					}	
				}
			}
	
}


	//*********************************************************************************************


	/// Fluid -- fill shell
	/// Include the modes of the outer surface, and exclude that of inner surface.
void IncVF::Fill_shell_vorticity(int shell_from_index)						
{
	
	TinyVector<complx,3> vorticity;
	complx vort_y;
	
	*G1 = 0.0;
	*G2 = 0.0;
	*G3 = 0.0;
	
	DP inner_radius = (*shell_radius)(shell_from_index - 1);
	DP outer_radius = (*shell_radius)(shell_from_index);
	
	DP kkmag;
	
	for (int l1=0;  l1<local_N1; l1++)
		for (int l2=0;  l2<N[2]; l2++) 
			for (int l3=0;  l3<=N[3]/2; l3++) 
			{
				kkmag = Kmagnitude(basis_type, l1, l2, l3, N, kfactor);
				
				if ( (kkmag > inner_radius) && (kkmag <= outer_radius) ) 
				{
					if (N[2] > 1)
					{	
						Compute_Modal_vorticity(basis_type, l1, l2, l3, N, 
												*V1, *V2, *V3, kfactor, vorticity);
						
						(*G1)(l1, l2, l3) = vorticity(0);
						(*G2)(l1, l2, l3) = vorticity(1);
						(*G3)(l1, l2, l3) = vorticity(2);
					}	
					
					else if (N[2] ==1)
					{
						Compute_Modal_vorticity_y_component(basis_type, l1, 0, l3, N, 
															*V1, *V1, *V3, kfactor, vort_y);
						
						(*G1)(l2, 0, l3) = vort_y;						
					}	
				}
			}

}



//*****************************  End of fill_sphere_shell.cc **********************************




