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
void IncVF::Compute_flux()
{	
	
	(*flux_self) = 0.0;
	
	if (ET_real_imag_switch == 1)
	{
		(*flux_self_real) = 0.0;
		(*flux_self_imag) = 0.0;
	}
	
	DP tot_real, tot_imag;

	for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
	{	
		Fill_in_sphere(sphere_index);	
					
		EnergyTr_Compute_nlin();											
		// nlin = U.grad U<	
		
		(*flux_self)(sphere_index) = - Prod_out_sphere_nlinV(sphere_index);	
		// flux_self = -(U.grad U<). U>
		
		if (ET_real_imag_switch == 1) 
		{	
			Prod_out_sphere_nlinV_real_imag(sphere_index, tot_real, tot_imag);
			(*flux_self_real)(sphere_index) = -tot_real;
			(*flux_self_imag)(sphere_index) = -tot_imag;	
		}						
	}
}


//*********************************************************************************************

// Scalar
//

void IncVF::Compute_flux(IncSF& T)
{

	(*flux_SF) = 0.0;
	
	if (ET_real_imag_switch == 1)
	{
		(*flux_SF_real) = 0.0;
		(*flux_SF_imag) = 0.0;
	}
	
	
	DP tot_real, tot_imag;
	
	Compute_flux();									
	// flux_self = (U.grad U<). U>	
	
	for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
	{		
		Fill_in_sphere(sphere_index, T);	
					
		EnergyTr_Compute_nlin(T);												
		// nlin = U.grad T<	
		
		(*flux_SF)(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, T);		
		// T.flux = (U.grad T<). T>		
		
		if (ET_real_imag_switch == 1) 
		{	
			Prod_out_sphere_nlinV_real_imag(sphere_index, T, tot_real, tot_imag);
			(*flux_SF_real)(sphere_index) = -tot_real;
			(*flux_SF_imag)(sphere_index) = -tot_imag;	
		}						
	}
}


//*********************************************************************************************

//
// Vector
void IncVF::Compute_flux(IncVF& W)
{
	
	(*flux_VF_in_out) = 0.0;
	(*flux_VF_in_in) = 0.0;
	
	if (ET_real_imag_switch == 1)
	{
		(*flux_VF_in_out_real) = 0.0;
		(*flux_VF_in_out_imag) = 0.0;
		
		(*flux_VF_in_in_real) = 0.0;
		(*flux_VF_in_in_imag) = 0.0;
	}
	
	
	DP tot_real, tot_imag;

	Compute_flux();									
	// flux_self = (U.grad U<). U>
	
	for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
	{		
		Fill_in_sphere(sphere_index);				
		
		EnergyTr_Compute_nlin(W);												
		// nlin = W.grad V<
		
		(*flux_VF_in_out)(sphere_index) = Prod_out_sphere_nlinV(sphere_index, W);	
		// (W.graad U<). W>
		
		(*flux_VF_in_in)(sphere_index) = Prod_in_sphere_nlinV(sphere_index, W);		
		// (W.graad U<). W<
		
		if (ET_real_imag_switch == 1) 
		{	
			Prod_out_sphere_nlinV_real_imag(sphere_index, W, tot_real, tot_imag);		
			(*flux_VF_in_out_real)(sphere_index) = tot_real;
			(*flux_VF_in_out_imag)(sphere_index) = tot_imag;
		
			Prod_in_sphere_nlinV_real_imag(sphere_index, W, tot_real, tot_imag);
			(*flux_VF_in_in_real)(sphere_index) = tot_real;
			(*flux_VF_in_in_imag)(sphere_index) = tot_imag;
		}	
	}
	
	
	
	
	// W< to W>
	
	(*W.flux_self) = 0.0;
	
	if (ET_real_imag_switch == 1)
	{
		(*W.flux_self_real) = 0.0;
		(*W.flux_self_imag)= 0.0;
	}
	
	for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
	{	
		Fill_in_sphere(sphere_index, W);										
		// G = W<
		
		EnergyTr_Compute_nlin();												
		// nlin = U.grad W<
		
		(*W.flux_self)(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, W);		
		// -(U.graad W<). W>
		
		if (ET_real_imag_switch == 1) 
		{	
			Prod_out_sphere_nlinV_real_imag(sphere_index, W, tot_real, tot_imag);		
			(*W.flux_self_real)(sphere_index) = -tot_real;
			(*W.flux_self_imag)(sphere_index) = -tot_imag;	
		}		
	}
	
	
	
	// W< to U>
	
	(*W.flux_VF_in_out) = 0.0;
	(*W.flux_VF_in_in) = 0.0;
	
	if (ET_real_imag_switch == 1)
	{
		(*W.flux_VF_in_out_real) = 0.0;
		(*W.flux_VF_in_out_imag) = 0.0;
		
		(*W.flux_VF_in_in_real) = 0.0;
		(*W.flux_VF_in_in_imag) = 0.0;
	}
	
	
	for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
	{	
		Fill_in_sphere(sphere_index, W);										
		// G = W<
		
		EnergyTr_Compute_nlin();												
		// nlin = U.grad W<
		
		(*W.flux_VF_in_out)(sphere_index) = Prod_out_sphere_nlinV(sphere_index);	
		//  (U.graad W<). U>
		
		(*W.flux_VF_in_in)(sphere_index) = Prod_in_sphere_nlinV(sphere_index);		
		//  (U.graad W<). U<
		
		if (ET_real_imag_switch == 1) 
		{	
			Prod_out_sphere_nlinV_real_imag(sphere_index, tot_real, tot_imag);		
			(*W.flux_VF_in_out_real)(sphere_index) = tot_real;
			(*W.flux_VF_in_out_imag)(sphere_index) = tot_imag;
		
			Prod_in_sphere_nlinV_real_imag(sphere_index, tot_real, tot_imag);
			(*W.flux_VF_in_in_real)(sphere_index) = tot_real;
			(*W.flux_VF_in_in_imag)(sphere_index) = tot_imag;
		}			
	}
	
	
	// U> to W>
	
	(*flux_VF_out_out) = 0.0;
	
	if (ET_real_imag_switch == 1)
	{
		(*flux_VF_out_out_real) = 0.0;
		(*flux_VF_out_out_imag) = 0.0;
	}
	
	
	for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
	{		
		Fill_out_sphere(sphere_index);										
		// G = U>
		
		EnergyTr_Compute_nlin(W);											
		// nlin = W.grad U>
		
		(*flux_VF_out_out)(sphere_index) = Prod_out_sphere_nlinV(sphere_index, W);	
		// (W.graad U>). W>
		
		if (ET_real_imag_switch == 1) 
		{	
			Prod_out_sphere_nlinV_real_imag(sphere_index, W, tot_real, tot_imag);
			(*flux_VF_out_out_real)(sphere_index) = tot_real;
			(*flux_VF_out_out_imag)(sphere_index) = tot_imag;
		}	
	}
	
	// Flux for Elsasser vars
	//
	
	(*flux_Elsasser) = 0.0;
	(*W.flux_Elsasser) = 0.0;
	
	if (ET_real_imag_switch == 1)
	{
		(*flux_Elsasser_real)= 0.0;
		(*flux_Elsasser_imag)= 0.0;
		
		(*W.flux_Elsasser_real)= 0.0;
		(*W.flux_Elsasser_imag)= 0.0;
	}
	
	
	IncVF::UB_to_Elsasser_field(W);												
	// U=Zp=(U+B); B=Zm=(U-B);
	
	// Flux: Zp to Zp
	for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
	{		
		Fill_in_sphere(sphere_index);											
		// G = Zp<
		
		EnergyTr_Compute_nlin(W);												
		// nlin = Zm.grad Zp<
		
		(*flux_Elsasser)(sphere_index) = -Prod_out_sphere_nlinV(sphere_index);	
		// (Zm.graad Zp<). Zp>
		
		if (ET_real_imag_switch == 1) 
		{	
			Prod_out_sphere_nlinV_real_imag(sphere_index, tot_real, tot_imag);
			(*flux_Elsasser_real)(sphere_index) = -tot_real;
			(*flux_Elsasser_imag)(sphere_index) = -tot_imag;
		}	
	}	
	
	
	// Flux: Zm to Zm
	for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
	{
		Fill_in_sphere(sphere_index, W);										
		// G = Zm<
		
		EnergyTr_Compute_nlin();												
		// nlin = Zp.grad Zm<
		
		(*W.flux_Elsasser)(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, W);	
		// (Zp.graad Zm<). Zm>
		
		if (ET_real_imag_switch == 1) 
		{	
			Prod_out_sphere_nlinV_real_imag(sphere_index, W, tot_real, tot_imag);
			(*W.flux_Elsasser_real)(sphere_index) = -tot_real;
			(*W.flux_Elsasser_imag)(sphere_index) = -tot_imag;
		}	
	}																																			// Flux: Zp to Zp
	
	IncVF::Elsasser_to_UB_field(W);													
	// Back to U, B vars
}


//*********************************************************************************************
//
// Magnetoconvection

void IncVF::Compute_flux(IncVF& W, IncSF& T)
{
	
	(*flux_SF) = 0.0;
	
	if (ET_real_imag_switch == 1)
	{
		(*flux_SF_real) = 0.0;
		(*flux_SF_imag) = 0.0;
	}
	
	
	DP tot_real, tot_imag;
	
	Compute_flux(W);
	
	for (int sphere_index = 1; sphere_index <= no_spheres; sphere_index++) 
	{		
		Fill_in_sphere(sphere_index, T);	
					
		EnergyTr_Compute_nlin(T);											
		// nlin = U.grad T<	
		
		(*flux_SF)(sphere_index) = -Prod_out_sphere_nlinV(sphere_index, T);		
		// T.flux = (U.grad T<). T>		
		
		if (ET_real_imag_switch == 1) 
		{			
			Prod_out_sphere_nlinV_real_imag(sphere_index, T, tot_real, tot_imag);
			(*flux_SF_real)(sphere_index) = -tot_real;
			(*flux_SF_imag)(sphere_index) = -tot_imag;
		}	
	}
}


//*********************************************************************************************
//
// RB convection

void IncVF::Compute_flux(IncSF& T, string Pr_switch)
{
	if (Pr_switch == "PRZERO") 
		Compute_flux();
	else
		Compute_flux(T);	
}

//
//

void IncVF::Compute_flux(IncVF& W, IncSF& T, string Pr_switch)
{
	if (Pr_switch == "PRZERO") 
		Compute_flux(W);
	else
		Compute_flux(W, T);
}		

//****************************  End of compute_flux.cc ****************************************






