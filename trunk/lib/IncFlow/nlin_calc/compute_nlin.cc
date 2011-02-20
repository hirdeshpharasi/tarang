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


/*! \file  compute_nlin.cc
 * 
 * @brief  Compute nonlinear \f$ N_i = \mathcal{F} (D_j V_j V_i) \f$.
 *
 *	Steps <BR>
 * # Inverse transform of  V: \f$ \vec{Vr} = \mathcal{F}(\vec{V}) \f$  <BR>
 * # nlin \f$ N_i  = (Vr_i)^2 \f$  <BR>
 * # nlin \f$ N_i  = D_i \mathcal{F}((Vr_i)^2) \f$ <BR>
 * # compute off-diagonal terms and put them in \f$ \vec{Vr} \f$. <BR>
 * # Forward transform \f$ \vec{Vr} \f$. <BR>
 * # Multiply by appropriate component of wavenumber to obtain 
 *		\f$ N_i = \mathcal{F} (D_j V_j V_i) \f$. The new terms are added to \f$ \vec{N} \f$
 *		that already contains diagonal terms. <BR>
 * # Dealias \f$ \vec{N} \f$ if DEALIASE switch is on. <BR>
 *
 *  @sa RSprod.cc
 *  @sa compute_derivative.cc
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug   TRANSPOSE order needs to be tested.
 */

#include "../IncVF.h"
#include "../IncSF.h"

//*********************************************************************************************
	    
void IncVF::Compute_nlin() 
{

// Normal order
//

#ifndef TRANSPOSE					
	
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3;	

	RV_Inverse_transform(*VF_temp_r);	
	
	// Vr[i] -> Vr[i]^2 stored in nlin[i]				                      
	Compute_RSprod_diag();										
	 
	// nlin[i]= Di F[Vr[i]^2] 
	NLIN_diag_Forward_transform_derivative(*VF_temp_r);			
  
  // Vr[i] = Vr[i]*Vr[j]; rules in RSprod
	Compute_RSprod_offdiag();	
	
	// Vr[i] = T(Vr[i]*Vr[j])																
	RV_Forward_transform_RSprod(*VF_temp_r);					

	// nlin[i] = Dj T(Uj * Ui)]
	Derivative_RSprod_VV();		
	
//
//		Transpose order
//

#else		

	 // *vir = Inv_tr(*Vi) 
	RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp);                     
	
	// nlin[i]= Di T[Vr[i]^2]
	Compute_RSprod_diag_ft_derivative();					
 
	// Vr[i] = Vr[i]*Vr[j]; rules in RSprod
	Compute_RSprod_offdiag();									
	
	// nlin[i] = Dj T(Uj * Ui)]	
	Forward_transform_derivative_RSprod_offdiag();				
	
#endif	

	if (alias_switch == "DEALIAS")		Dealias_nlin();
}   


//*********************************  End of compute_nlin.cc  **********************************




