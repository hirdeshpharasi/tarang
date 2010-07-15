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
/*! \file  prod_nlinV.cc
 * 
 * @brief  Returns \f$\sum \Re(N_i^*(k) V_i(k)) \f$  or \f$ \sum \Re(N_i^*(k) W_i(k)) \f$
 *			for the specified region  (inside sphere or outside sphere).  
 *			\f$ \vec{N} \f$ is the nlin term. 
 * 
 * We also compute ans_real = \f$ \sum \Re(N_i(k)) * \Re(V_i(k)) \f$ and 
 *				   ans_imag = \f$ \sun \Im(N_i(k)) * \Im(V_i(k)) \f$ 
 *					(Could be \fR W_i \f$ as well).
 *
 * Sphere_index is passed to the function as a parameter.  The function
 *		picks the radius and computes real( A[k]*conj(B[k]) ) using shell-mult functions. <BR>
 *		For inside sphere, the surface modes are included.
 *		Fou outside sphere, the surface modes are excluded.
 *
 * @sa void Shell_mult_single(...)
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

/** @brief For fluid: compute \f$\sum \Re(N_i^*(k) V_i(k)) \f$ outside sphere.  
 *	
 *	@param	sphere_index 1:no_sphere, the last sphere contains all the modes.
 *
 *	@return  \f$\sum \Re(N_i^*(k) V_i(k)) \f$ outside sphere.
 */
DP IncVF::Prod_out_sphere_nlinV(int sphere_index)
{
	DP	radius = (*sphere_radius)(sphere_index);			// radius of the sphere	

	return  ( Shell_mult_single(basis_type, alias_switch, N, *nlin1, *V1, 
								radius, INF_RADIUS, kfactor)
			+ Shell_mult_single(basis_type, alias_switch, N, *nlin2, *V2, 
								radius, INF_RADIUS, kfactor)
			+ Shell_mult_single(basis_type, alias_switch, N, *nlin3, *V3, 
								radius, INF_RADIUS, kfactor)); 
}

// FOR PASSIVE SCALAR --- OUTSIDE SPHERE

DP IncVF::Prod_out_sphere_nlinV(int sphere_index, IncSF& T)
{
	DP	radius = (*sphere_radius)(sphere_index);			
	
	return ( Shell_mult_single(basis_type, alias_switch, N, *nlin1, *T.F, 
								radius, INF_RADIUS, kfactor) );			
}


// FOR MHD 

DP IncVF::Prod_out_sphere_nlinV(int sphere_index, IncVF& W)
{
	DP	radius = (*sphere_radius)(sphere_index);				

	return  ( Shell_mult_single(basis_type, alias_switch, N, *nlin1, *W.V1, 
								radius, INF_RADIUS, kfactor)
			+ Shell_mult_single(basis_type, alias_switch, N, *nlin2, *W.V2, 
								radius, INF_RADIUS, kfactor)
			+ Shell_mult_single(basis_type, alias_switch, N, *nlin3, *W.V3, 
								radius, INF_RADIUS, kfactor));
}


//*********************************************************************************************

/** @brief For fluid: compute \f$\sum \Re(N_i^*(k) V_i(k)) \f$ inside sphere.  
 *	
 *	@param	sphere_index 1:no_sphere, the last sphere contains all the modes.
 *
 *	@return  \f$\sum \Re(N_i^*(k) V_i(k)) \f$ inside sphere.
 */
DP IncVF::Prod_in_sphere_nlinV(int sphere_index)
{
	DP	radius = (*sphere_radius)(sphere_index);			// radius of the sphere	

	return  ( Shell_mult_single(basis_type, alias_switch, N, *nlin1, *V1, 0, radius, kfactor)
			+ Shell_mult_single(basis_type, alias_switch, N, *nlin2, *V2, 0, radius, kfactor)
			+ Shell_mult_single(basis_type, alias_switch, N, *nlin3, *V3, 0, radius, kfactor)); 
}

// FOR PASSIVE SCALAR --- INSIDE SPHERE

DP IncVF::Prod_in_sphere_nlinV(int sphere_index, IncSF& T)
{
	DP	radius = (*sphere_radius)(sphere_index);			// radius of the sphere	
	
	return (Shell_mult_single(basis_type, alias_switch, N, *nlin1, *T.F, 0, radius, kfactor) );	
}


// FOR MHD 

DP IncVF::Prod_in_sphere_nlinV(int sphere_index, IncVF& W)
{
	DP	radius = (*sphere_radius)(sphere_index);			// radius of the sphere	

	return(Shell_mult_single(basis_type, alias_switch, N, *nlin1, *W.V1, 0, radius, kfactor)
		 + Shell_mult_single(basis_type, alias_switch, N, *nlin2, *W.V2, 0, radius, kfactor)
		 + Shell_mult_single(basis_type, alias_switch, N, *nlin3, *W.V3, 0, radius, kfactor));
}


//*********************************************************************************************

//
// Real*real & imag*imag parts computed separately.
//

// FOR FLUID --- OUTSIDE SPHERE

void IncVF::Prod_out_sphere_nlinV_real_imag(int sphere_index, DP& tot_real, DP& tot_imag)
{
	DP	radius = (*sphere_radius)(sphere_index);			// radius of the sphere	
	DP  temp_real = 0.0;
	DP	temp_imag = 0.0;
	
	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin1, *V1, 
								radius, INF_RADIUS, tot_real, tot_imag, kfactor);
								
	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin2, *V2, 
								radius, INF_RADIUS, temp_real, temp_imag, kfactor);
								
	tot_real += temp_real;		
	tot_imag += temp_imag;
	
	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin3, *V3, 
								radius, INF_RADIUS, temp_real, temp_imag, kfactor);
								
	tot_real += temp_real;		
	tot_imag += temp_imag;
}

// Scalar

void IncVF::Prod_out_sphere_nlinV_real_imag
(
	int sphere_index, 
	IncSF& T, 
	DP& tot_real, DP& tot_imag
)
{
	DP	radius = (*sphere_radius)(sphere_index);			// radius of the sphere	
	
	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin1, *T.F, 
								radius, INF_RADIUS, tot_real, tot_imag, kfactor);
}

// MHD
void IncVF::Prod_out_sphere_nlinV_real_imag
(
	int sphere_index, 
	IncVF& W, 
	DP& tot_real, DP& tot_imag
)
{
	DP	radius = (*sphere_radius)(sphere_index);			// radius of the sphere	
	DP  temp_real = 0.0;
	DP	temp_imag = 0.0;
	
	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin1, *W.V1, 
								radius, INF_RADIUS, tot_real, tot_imag, kfactor);
								
	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin2, *W.V2, 
								radius, INF_RADIUS, temp_real, temp_imag, kfactor);
								
	tot_real += temp_real;		
	tot_imag += temp_imag;
	
	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin3, *W.V3, 
								radius, INF_RADIUS, temp_real, temp_imag, kfactor);
								
	tot_real += temp_real;		
	tot_imag += temp_imag;
}


//*********************************************************************************************

//
// Real*real & imag*imag parts computed separately.
//

// FOR FLUID --- INSIDE SPHERE

void IncVF::Prod_in_sphere_nlinV_real_imag(int sphere_index, DP& tot_real, DP& tot_imag)
{
	DP	radius = (*sphere_radius)(sphere_index);			// radius of the sphere	
	DP  temp_real = 0.0;
	DP	temp_imag = 0.0;
	

	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin1, *V1, 0, 
								radius, tot_real, tot_imag, kfactor);
								
	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin2, *V2, 0, 
								radius, temp_real, temp_imag, kfactor);
								
	tot_real += temp_real;		
	tot_imag += temp_imag;
	
	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin3, *V3, 0, 
								radius, temp_real, temp_imag, kfactor);
								
	tot_real += temp_real;		
	tot_imag += temp_imag;
}

// Scalar

void IncVF::Prod_in_sphere_nlinV_real_imag
(
	int sphere_index, 
	IncSF& T, 
	DP& tot_real, DP& tot_imag
)
{
	DP	radius = (*sphere_radius)(sphere_index);			// radius of the sphere	
	
	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin1, *T.F, 0, 
								radius, tot_real, tot_imag, kfactor);
}

// MHD
void IncVF::Prod_in_sphere_nlinV_real_imag
(
	int sphere_index, 
	IncVF& W, 
	DP& tot_real, DP& tot_imag
)
{
	DP	radius = (*sphere_radius)(sphere_index);			// radius of the sphere	
	DP  temp_real = 0.0;
	DP	temp_imag = 0.0;

	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin1, *W.V1, 0, 
								radius, tot_real, tot_imag, kfactor);
								
	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin2, *W.V2, 0, 
								radius, temp_real, temp_imag, kfactor);
								
	tot_real += temp_real;		
	tot_imag += temp_imag;
	
	Shell_mult_single_real_imag(basis_type, alias_switch, N, *nlin3, *W.V3, 0, 
								radius, temp_real, temp_imag, kfactor);
								
	tot_real += temp_real;		
	tot_imag += temp_imag;
}


//*********************************************************************************************
// Helicity flux
DP IncVF::Prod_out_sphere_nlin_vorticity(int sphere_index)
{
	DP	radius = (*sphere_radius)(sphere_index);			// radius of the sphere	
	
	return Shell_mult_vorticity(basis_type, alias_switch, N, 
								*nlin1, *nlin2, *nlin3, *V1, *V2, *V3,
								radius, INF_RADIUS, kfactor);
	
}


	// Helicity flux
DP IncVF::Prod_out_sphere_nlin_vector_potential(int sphere_index, IncVF& W)
{
	DP	radius = (*sphere_radius)(sphere_index);			// radius of the sphere	
	
	return Shell_mult_vector_potential(basis_type, alias_switch, N, 
									   *nlin1, *nlin2, *nlin3, *W.V1, *W.V2, *W.V3,
									   radius, INF_RADIUS, kfactor);
	
}

//******************************  End of Prod_nlinV.cc ****************************************






