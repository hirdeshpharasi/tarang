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


/*! \file  compute_nlin_VF_SF.cc
 * 
 * @brief  Compute nonlinear \f$ N_i^V = \mathcal{F} (D_j (V_j V_i - B_j B_i)) \f$ and 
 *		\f$ N_i^B = \mathcal{F} (D_j (V_j B_i - B_j V_i)) \f$,
 *		\f$ N_i^\psi = \mathcal{F} (D_j V_j \psi) \f$.
 *
 *  We use compute_nlin(W) to compute MHD nlin terms, and use the method of 
 *	compute_nlin(T) to compute the temperature nlin term.
 *
 *	@sa compute_nlin_VF_MHD.cc
 *  @sa compute_nlin_SF.cc
 *  @sa RSprod.cc
 *  @sa compute_derivative.cc
 *
 * @note  The basis functions of V and B are the same in the present implementation.
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Sept 2008
 *
 * @bug  Transpose order
 */


#include "../IncVF.h"
#include "../IncSF.h"


//*********************************************************************************************

void IncVF::Compute_nlin(IncVF& W, IncSF& T) 
{

	// compute U.nlin and W(B).nlin
	IncVF:: Compute_nlin(W);							

//
//	Now compute T.nlin = Dj T [Vj T]	
//

#ifndef TRANSPOSE

	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3; 
	
	*T.Fr = *T.F;			
	
	// Vr[i]= I(Vk[i])
	RV_Inverse_transform(*VF_temp_r);							
  
	// Compute T.nlin = Dj FT(Vj*T) 
	
	// Inv Transform of *T.Fr 
	T.RS_Inverse_transform(*VF_temp_r);								

	// T.nlin = V1r * (T.Fr)
	Array_real_mult(N, *V1r, *T.Fr, *T.nlin);				
	
	// T.SF_temp = V2r * (T.Fr)
	Array_real_mult(N, *V2r, *T.Fr, *T.SF_temp);			
		
	Forward_transform_array(basis_type, N, *T.nlin, *VF_temp_r, 0);	 
	// parity = 0; sin*sin =cos
	
	Forward_transform_array(basis_type, N, *T.SF_temp, *VF_temp_r, 1);	
	// parity = 1; sin*cos =sin		

	Xderiv_RSprod_VT(*T.nlin, *T.nlin);					 
	Yderiv_RSprod_VT(*T.SF_temp, *T.SF_temp);					 
    
	// T.nlin= Dj T[Vj F]     
	*T.nlin = *T.nlin + *T.SF_temp;							

	Array_real_mult(N, *V3r, *T.Fr, *T.SF_temp); 
	
	Forward_transform_array(basis_type, N, *T.SF_temp, *VF_temp_r, 1);	
	// parity = 1; sin*cos =sin
	
	Zderiv_RSprod_VT(*T.SF_temp, *T.SF_temp);
  
	*T.nlin = *T.nlin + *T.SF_temp; 
		
//*********************************************************************************************
//  TRANSPOSE_ORDER 
	
#else		

	RVF::RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp);   
	// *Vir = Inv_tr(*Vi) 
	         
	T.RS_Inverse_transform_transpose_order(*F, *VF_temp);							
	// *VF_temp = Inv_tr(*F) 
	
	Array_real_mult(N, *V1r, *T.Fr, *VF_temp_r);
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *T.nlin, 0);
	Xderiv_RSprod_VT(*T.nlin, *T.nlin);	
									
	Array_real_mult(N, *V2r, *T.Fr, *VF_temp_r);			
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *VF_temp, 1);
	Yderiv_RSprod_VT(*VF_temp, *VF_temp);	
	*T.nlin = *T.nlin + *VF_temp;
	
	Array_real_mult(N, *V3r, *T.Fr, *VF_temp_r);			
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *VF_temp, 1);
	Zderiv_RSprod_VT(*VF_temp, *VF_temp);	
	*T.nlin = *T.nlin + *VF_temp;

#endif		

	if (alias_switch == "DEALIAS")		Dealias_nlin();
	
}

//*********************************************************************************************
// For RB convection
//

void IncVF::Compute_nlin(IncVF& W, IncSF& T, string Pr_switch) 
{

	if (Pr_switch == "PRZERO") 
		Compute_nlin(W);
	else
		Compute_nlin(W, T);
}		


//***************************** fn compute_nlin_VF_SF ends ************************************


