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

/*! \file  compute_nlin_SF.cc
 * 
 * @brief  Compute nonlinear \f$ N_i^V = \mathcal{F} (D_j V_j V_i) \f$ and 
 *		\f$ N_i^\psi = \mathcal{F} (D_j V_j \psi) \f$.
 *
 *	Steps <BR>
 * # Inverse transform of V & T: \f$ \vec{Vr} = \mathcal{F}(\vec{V}) \f$  & 
 *			\f$ Fr = \mathcal{F}(F) \f$  <BR>
 * # T.nlin \f$ N^\psi = \mathcal{F} (D_j V_j \psi) \f$. <BR>
 * # nlin \f$ N_i  = (Vr_i)^2 \f$  <BR>
 * # nlin \f$ N_i  = \mathcal{F}(D_i (Vr_i)^2) \f$ <BR>
 * # compute off-diagonal terms and put them in \f$ \vec{Vr} \f$. <BR>
 * # Forward transform \f$ \vec{Vr} \f$. <BR>
 * # Multiply by appropriate component of wavenumber to obtain 
 *		\f$ N_i = \mathcal{F} (D_j V_j V_i) \f$. The new terms are added to \f$ \vec{N} \f$
 *		that already contains diagonal terms. <BR>
 * # Dealias \f$ \vec{N} \f$ if DEALIASE switch is on. <BR>
 *  
 * # RB Convection: If Pr==0: nlin(), else nlin(T).
 *
 *  @sa RSprod.cc
 *  @sa compute_derivative.cc
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug   Transpose order needs to tested
 */

  
#include "../IncVF.h"
#include "../IncSF.h"



//**************************************************************************************************


void IncVF::Compute_nlin_NonBoussinesq(IncSF& T) 
{
	
	//	Normal order
	//
	
	
#ifndef TRANSPOSE
	
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3;	
	
	*T.Fr = *T.F;	
	
	RVF::RV_Inverse_transform(*VF_temp_r);							
	
	// Compute T.nlin = Dj FT(Vj*T) 
	
	T.RS_Inverse_transform(*VF_temp_r);								// Inv Transform of *T.Fr 
	
	Array_real_mult(N, *V1r, *T.Fr, *T.nlin);						// T.nlin = V1r * (T.Fr)
	Array_real_mult(N, *V2r, *T.Fr, *T.SF_temp);					// T.SF_temp = V2r * (T.Fr)
	
	Forward_transform_array(basis_type, N, *T.nlin, *VF_temp_r, 0);	    
	// parity = 0; sin*sin =cos
	
	Forward_transform_array(basis_type, N, *T.SF_temp, *VF_temp_r, 1);	
	// parity = 1; sin*cos =sin	
	
	
	Xderiv_RSprod_VT(*T.nlin, *T.nlin);					 
	Yderiv_RSprod_VT(*T.SF_temp, *T.SF_temp);	
	
	
	*T.nlin = *T.nlin + *T.SF_temp;									// T.nlin= Dj T[Vj  F]
	
	Array_real_mult(N, *V3r, *T.Fr, *T.SF_temp); 
	
	Forward_transform_array(basis_type, N, *T.SF_temp, *VF_temp_r, 1);	
	// parity = 1; sin*cos =sin
	
	Zderiv_RSprod_VT(*T.SF_temp, *T.SF_temp);
	
	*T.nlin = *T.nlin + *T.SF_temp; 
	
	
	// End of Scalar field computation 
	
	// Now Compute Dj (totrho Ui Uj)
	
	// T.Compute_total_density()::: T.Fr = rho0+(1-x)+rho'
	//T.Fr contains density in real space
	// total_rho = rho'+mean density + rho(y); rho(x) = 1-x
	int kx;
	DP x;
	for (int l1=0; l1<local_N1; l1++) {
		kx = Get_kx("SCFT", l1, N);
		x = (kx+0.5)/N[1];
		
		for (int l2=0; l2<N[2]; l2++)
			for (int l3=0; l3<(N[3]/2); l3++) 
				(*T.Fr)(l1, l2, l3) += (-x);
	}
	
	real(*T.Fr) = globalvar_alpha_DT*real(*T.Fr) + 1.0; 
	imag(*T.Fr) = globalvar_alpha_DT*imag(*T.Fr) + 1.0;  
	// mean density = 1.0
	// *T.Fr contains total rho/rho_0 now. 

	//	Compute_RSprod_diag():   Vr[i] -> rho Vr[i]^2 stored in nlin[i)
	
	Array_real_mult(N, *V1r, *V1r, *nlin1);  
	Array_real_mult(N, *V2r, *V2r, *nlin2);
	Array_real_mult(N, *V3r, *V3r, *nlin3);
	
	Array_real_mult(N, *nlin1, *T.Fr, *nlin1);  
	Array_real_mult(N, *nlin2, *T.Fr, *nlin2);
	Array_real_mult(N, *nlin3, *T.Fr, *nlin3); 
	
	
	NLIN_diag_Forward_transform_derivative(*VF_temp_r);	
	// nlin[i]=  Di T[rho Vr[i]^2.. we assume that rho0 >> rho' so the parity is still 
	// SIN for nlin1 and COS for other 2.
	
	// Compute_RSprod_offdiag() :::		Vr[i]=rho * Vr[i]*Vr[j];
	*VF_temp = *V1r;  
	
	Array_real_mult(N, *V1r, *V2r, *V1r);		// V1r= V1r*V2r
	Array_real_mult(N, *V2r, *V3r, *V2r);		// V2r= V2r*V3r
	Array_real_mult(N, *V3r, *VF_temp, *V3r);	// V3r= V3r*V1r 
	
	Array_real_mult(N, *V1r, *T.Fr, *V1r);  
	Array_real_mult(N, *V2r, *T.Fr, *V2r);
	Array_real_mult(N, *V3r, *T.Fr, *V3r);
	
	RV_Forward_transform_RSprod(*VF_temp_r);			// Vr[i]=T(rho Vr[i]*Vr[j])
	
	Derivative_RSprod_VV();								// nlin[i] = Dj T( Uj Ui)
	
	//*********************************************************************************************	
	// TRANSPOSE ORDER
	
#else		
	
	RVF::RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp);  
	// *Vir = Inv_tr(*Vi) 
	
	T.RS_Inverse_transform_transpose_order(*F, *VF_temp);							
	// *T.Fr = Inv_tr(*F) ; *VF_temp is temporary storage
	
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
	
	// End of Scalar field computation 
	
	// Now for the velocity field
	
	// total_rho = rho'+ rho(x) + mean density(=1); rho(x) = -x 
	// delta_rho/rho = alpha DT * (T_cond + theta)
	// Density at the bottom plate is rho_0 (mean density). Here normalized to 1.
	
	int kx;
	DP x;
	for (int l1=0; l1<local_N1; l1++) {
		kx = Get_kx("SCFT", l1, N);
		x = (kx+0.5)/N[1];
		
		for (int l2=0; l2<N[2]; l2++)
			for (int l3=0; l3<(N[3]/2); l3++)
				(*T.Fr)(l1, l2, l3) += (-x);
	}
	
	*T.Fr = globalvar_alpha_DT*(*T.Fr) + 1.0;  // mean density = 1.0
	
	// *T.Fr contains total rho/rho_0 now. 
	
	// Compute_RSprod_diag_ft_derivative();	nlin[i]= Di T[rho Vr[i]^2]
	
	Array_real_mult(N, *V1r, *V1r, *VF_temp_r);
	Array_real_mult(N, *VF_temp_r, *T.Fr, *VF_temp_r);
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin1, 0);
	Xderiv(basis_type, N, *nlin1, *nlin1, kfactor, 0);
	
	Array_real_mult(N, *V2r, *V2r, *VF_temp_r);
	Array_real_mult(N, *VF_temp_r, *T.Fr, *VF_temp_r);
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin2, 0);
	Yderiv(basis_type, N, *nlin2, *nlin2, kfactor, 0);  // parity unimportant here.
	
	Array_real_mult(N, *V3r, *V3r, *VF_temp_r); 
	Array_real_mult(N, *VF_temp_r, *T.Fr, *VF_temp_r);
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin3, 0);
	Zderiv(basis_type, N, *nlin3, *nlin3, kfactor, 0);	 // parity unimportant here.
	
	// Vr[i] = rho Vr[i]*Vr[j]; rules in RSprod
	// IncVF::Compute_RSprod_offdiag();	
	
	*VF_temp = *V1r;  
	
	Array_real_mult(N, *V1r, *V2r, *V1r);		// V1r= V1r*V2r
	Array_real_mult(N, *V2r, *V3r, *V2r);		// V2r= V2r*V3r
	Array_real_mult(N, *V3r, *VF_temp, *V3r);	// V3r= V3r*V1r 
	
	Array_real_mult(N, *V1r, *T.Fr, *V1r);  
	Array_real_mult(N, *V2r, *T.Fr, *V2r);
	Array_real_mult(N, *V3r, *T.Fr, *V3r);
	
	
	
	IncVF::Forward_transform_derivative_RSprod_offdiag();		// nlin[i] = Dj T(Uj * Ui)]	
	
#endif	
	
	if (alias_switch == "DEALIAS")		Dealias_nlin(T);
	
}
//***************************** fn compute_nlin_SF ends ***************************************


