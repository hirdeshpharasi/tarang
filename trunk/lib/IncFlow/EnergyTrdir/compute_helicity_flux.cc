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

/*! \file  compute_flux.cc
 * 
 * @brief  Computes energy flux from inside/outside a sphere of field V/W to inside/outsude
 *			of field V/W.
 *
 *	The Giver field is filled inside/outside the sphere (region A1).  
 *	We do the same for the receiver field (region A2). <BR>
 * 
 *	Definitions given in EnergyTr.h.  We compute the energy transfer for these regions.
 *
 * @sa EnergyTr.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "../IncVF.h"
#include "../IncSF.h"

//*********************************************************************************************

// Fluid
void IncVF::Compute_kinetic_helicity_flux()
{	
	
	(*flux_hk) = 0.0;

	for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
	{	
		Fill_in_sphere(sphere_index);	
		
		EnergyTr_Compute_nlin();										
		// nlin = U.grad U<	
		
		(*flux_hk)(sphere_index) =  - Prod_out_sphere_nlin_vorticity(sphere_index)/2;	
		// flux_hm = -(U.grad U<). omega>
		
		
		Fill_in_sphere_vorticity(sphere_index);
		// G = i k x u(k)
		
		EnergyTr_Compute_nlin();										
		// nlin = U.grad omega	
		
		(*flux_hk)(sphere_index) +=   -Prod_out_sphere_nlinV(sphere_index)/2;	
		// Add -(U.grad omega). U>
		
		Fill_in_sphere(sphere_index);	
		
		EnergyTr_Compute_nlin_vorticity_helper();										
		// nlin = omega.grad U<	
		
		(*flux_hk)(sphere_index) +=   Prod_out_sphere_nlinV(sphere_index)/2;	
		// Add -(omega.grad U<). U>
		
							
	}
}


//*********************************************************************************************

//
// Vector
void IncVF::Compute_magnetic_helicity_flux(IncVF& W)
{
	
	(*W.flux_hk) = 0.0;			// magnetic helicity
	
	
	for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
	{
		
		Fill_in_sphere(sphere_index);								
		// G = U<
		
		EnergyTr_Compute_nlin(W);	
		// nlin = W.grad U<
		
		
		(*W.flux_hk)(sphere_index) = Prod_out_sphere_nlin_vector_potential(sphere_index, W)/2;	
		// flux_hm = W.grad U<. a>
		
		
		
		Fill_in_sphere(sphere_index, W);								
		// G = W<
		
		EnergyTr_Compute_nlin();	
		// nlin = U.grad W<
		
		(*W.flux_hk)(sphere_index) += -Prod_out_sphere_nlin_vector_potential(sphere_index, W)/2;	
		// add -(U.graad W<). a>
		
		
		
		Fill_in_sphere(sphere_index, W);								
		// G = W<
		
		EnergyTr_Compute_nlin_UcrossB();
		// nlin = u x W<
		
		(*W.flux_hk)(sphere_index) += Prod_out_sphere_nlinV(sphere_index, W)/2;		
		// add (U x W<). B>
		
	}
	
}



//*********************************************************************************************

// Fluid 2D
void IncVF::Compute_enstrophy_flux()
{	
	
	if (N[2] == 1)
	{
		(*flux_hk) = 0.0;
		
		
		for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
		{
			
			Fill_in_sphere_vorticity(sphere_index);
			// G = i k x u(k): 2D
			
			EnergyTr_Compute_nlin();										
			// nlin = U.grad omega	
			
			(*flux_hk)(sphere_index) =   - Prod_out_sphere_nlin_vorticity(sphere_index);	
			// -(U.grad omega). omega> : 2D
		}
	}
}


//*********************************************************************************************

// MHD 2D
void IncVF::Compute_magnetic_enstrophy_flux(IncVF& W)
{	
	
	if (N[2] == 1)
	{
		(*W.flux_hk) = 0.0;
		
		
		for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
		{
			
			Fill_in_sphere(sphere_index, W);
			// G = W<
			
			EnergyTr_Compute_nlin_UcrossB();
			// nlin = u x W< : 2D
			
			(*flux_hk)(sphere_index) = Prod_out_sphere_nlin_vector_potential(sphere_index, W);	
			// (U.grad omega). omega>
		}
	}
}


//****************************  End of compute_flux.cc ****************************************






