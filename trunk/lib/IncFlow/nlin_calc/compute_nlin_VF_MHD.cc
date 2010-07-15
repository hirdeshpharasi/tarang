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


/*! \file  compute_nlin_VF_MHD.cc
 * 
 * @brief  Compute nonlinear \f$ N_i^V = \mathcal{F} (D_j (V_j V_i - B_j B_i)) \f$ and 
 *		\f$ N_i^B = \mathcal{F} (D_j (V_j B_i - B_j V_i)) \f$.
 *
 *  We use \f$ \vec{Z}^\pm \f$ variables that decreases the no of required FFTs.
 *
 * @note \f$ N_i^V = \mathcal{F} (D_j (V_j V_i - B_j B_i)) 
 *				= \mathcal{F} (D_j (Z_j^- Z_i^+ + Z^+_j Z^-_i))\f$
 *
 * @note \f$ N_i^B = \mathcal{F} (D_j (V_j B_i - B_j V_i))
 *				= \mathcal{F} (D_j (Z_j^- Z_i^+ - Z^+_j Z^-_i))\f$
 *
 *	Steps <BR>
 *  # (V,B) <- \f$ \vec{Z}^+,  \vec{Z}^- \f$. <BR>
 *  # Inverse transform of \f$ \vec{Z}^\pm \f$ and store in real arrays. <BR>
 * # nlin \f$ N_i  = (Vr_i * Wr_i) \f$  <BR>
 * # nlin \f$ N_i  = D_i \mathcal{F}(Vr_i Wr_i) \f$ <BR>
 * # compute off-diagonal terms \f$ (Vr_i Wr_j) \f$ and put them in 
 *			\f$ \vec{Vr},  \vec{Wr}\f$. <BR>
 * # Forward transform \f$ \vec{Vr} \f$. <BR>
 * # Multiply by appropriate component of wavenumber to obtain 
 *		\f$ N_i = \mathcal{F} (D_j V_j W_i) \f$. The new terms are added to \f$ \vec{N} \f$
 *		that already contains diagonal terms. <BR>
 * # Go back (V,B) vars from \f$ (Z^+, Z^-) \f$.
 * # Dealias \f$ \vec{N} \f$ if DEALIASE switch is on. <BR>
 *
 *  @sa RSprod.cc
 *  @sa compute_derivative.cc
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Sept 2008
 *
 * @bug   Transpose order needs to be tested.
 */


#include "../IncVF.h"
#include "../IncSF.h"


//*********************************************************************************************	

void IncVF::Compute_nlin(IncVF& W) 
{

#ifndef TRANSPOSE

	// U=Zp=(U+B); B=Zm=(U-B);
	IncVF::UB_to_Elsasser_field(W);									
	

	*V1r = *V1;		
	*V2r = *V2;		
	*V3r = *V3;	
					
	*W.V1r = *W.V1; 
	*W.V2r = *W.V2; 
	*W.V3r =*W.V3;	
	
	
	RVF::RV_Inverse_transform(*VF_temp_r); 
	
	W.RV_Inverse_transform(*VF_temp_r);	

	// nlin[i] =  Vi Wi
	Compute_RSprod_diag(W);									
	
	// nlin[i]=	 T[Vi Wi]
	NLIN_diag_Forward_transform_derivative(*VF_temp_r);				
	
 
	*W.nlin1 = *nlin1;  
	*W.nlin2 = *nlin2;	
	*W.nlin3 = *nlin3;		
 
	// cross terms: Vj Wi	
	Compute_RSprod_offdiag(W);								
  
	// Vr[i] = T(Vj Wi)
	RV_Forward_transform_RSprod(*VF_temp_r);								
	
	W.RV_Forward_transform_RSprod(*VF_temp_r);	

	// Compute derivatives and construct nlin
	Derivative_RSprod_VV(W);									

	// Convert to V-B vars

	// U.nlin=(Zp.nlin+Zm.nlin)/2; B.nlin=(Zp.nlin-Zm.nlin)/2;	
	Elsasser_to_UB_nlin(W);	
		
	// Back to U, B vars
	Elsasser_to_UB_field(W);									


//*********************************************************************************************
//	Transpose Order 


#else

	// U=Zp=(U+B); B=Zm=(U-B);
	IncVF::UB_to_Elsasser_field(W);									

	
	RVF::RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp);
	
	W.RV_Inverse_transform_transpose_order(*W.V1, *W.V2, *W.V3, *W.VF_temp);
		  
	// nlin[i]= Di T[Vr[i]*Wr[i]]	  
	IncVF::Compute_RSprod_diag_ft_derivative(W);					


	*W.nlin1 = *nlin1;  
	*W.nlin2 = *nlin2;	
	*W.nlin3 = *nlin3;		
 
	// cross terms: Vr[1]=Vr[1]*Wr[2];  W1r = W1r*V2r...  
	IncVF::Compute_RSprod_offdiag(W);								
	
	IncVF::Forward_transform_derivative_RSprod_offdiag(W);
  

	// Convert to V-B vars

	// U.nlin=(Zp.nlin+Zm.nlin)/2; B.nlin=(Zp.nlin-Zm.nlin)/2;
	IncVF::Elsasser_to_UB_nlin(W);	
			
	// Back to U, B vars
	IncVF::Elsasser_to_UB_field(W);										

#endif	

	if (alias_switch == "DEALIAS")		Dealias_nlin(W);
}

//***************************** fn compute_nlin_VF_MHD ends ***********************************



