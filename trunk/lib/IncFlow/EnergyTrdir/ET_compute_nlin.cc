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

/*! \file  ET_compute_nlin.cc
 * 
 * @brief Compute energytransfer(ET) (a) nlin \f$ N_i \leftarrow D_j \mathcal{F}(H_j * G_i) \f$
 *			where G is the giver field, where H is the helper field (V or W); <BR>
 *			(b) nlin \f$ N_i \leftarrow D_j \mathcal{F}(Vr_j * T^<) \f$.
 *	
 *	Steps for VF -> VF <BR>
 * # Inverse transform of  H (Helper, V or W). <BR>
 *
 * # Inverse transform G (Giver). <BR>
 *
 * # nlin \f$ N_i  = H_i * G_i \f$  <BR>
 *
 * # nlin \f$ N_i  = D_i \mathcal{F}(H_i * G_i) \f$ <BR>
 *
 * # compute off-diagonal terms \f$ H_j * G_i \f$ and put them in \f$ \vec{Vr} \f$. <BR>
 *
 * # Forward transform \f$ \vec{Vr} \f$. <BR>
 *
 * # Multiply by appropriate component of wavenumber to obtain 
 *		\f$ N_i = \mathcal{F} (D_j V_j G_i) \f$. The new terms are added to \f$ \vec{N} \f$
 *		that already contains diagonal terms. <BR>
 * 
 *  
 *	Steps for SF -> SF <BR>
 *	# Inverse transform of  V (Helper). <BR>
 *
 *	# Store \f$ T<  \f$ in \f$ G_1 \f$.
 *
 *	# nlin \f$ N_i \leftarrow \mathcal{F}(Vr_i * G_1) \f$ <BR>
 *
 *	# Apply derivative: \f$ N_1 \leftarrow D_i \mathcal{F}(Vr_i * G_1) \f$ <BR>
 *
 *  # RB Convection: If Pr==0: nlin(), else nlin(T).
 *
 *  @sa RSprod.cc
 *  @sa compute_derivative.cc
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug  No known bugs
 */

#include "../IncVF.h"
#include "../IncSF.h"					


//*********************************************************************************************

/*! @brief Compute energytransfer(ET) nlin \f$ N_i \leftarrow D_j \mathcal{F}(Vr_j * G_i) \f$
 *			G is the giver field, i.e., truncated V itself.
 *
 *	@sa ET_compute_nlin.cc
 *
 *  @return \f$ N_i \leftarrow D_j \mathcal{F}(Vr_j * G_i) \f$.
 */ 
void IncVF::EnergyTr_Compute_nlin() 
{  

#ifndef TRANSPOSE
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3;     	          

	RV_Inverse_transform(*VF_temp_r);  
	
	ET_Inverse_transform_G(*VF_temp_r);
											  
	Compute_RSprod_diagET();					
 
	NLIN_diag_Forward_transform_derivative(*VF_temp_r);                     
  
	Compute_RSprod_offdiagET();	
	
	RV_Forward_transform_RSprod(*VF_temp_r);			// Vr[i]=FFT(Vr[i]*Vr[j])		
	ET_Forward_transform_RSprod_G(*VF_temp_r);
	
	Derivative_RSprod_VV_ET();	
		

#else	
	RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp);
	
	ET_Inverse_transform_G_transpose_order();
						                      
	Compute_RSprod_diag_ft_derivativeET();			
 
	Compute_RSprod_offdiagET();							// Now Vir, Gi contains cross products
					
	Forward_transform_derivative_RSprod_offdiagET();		
#endif			
}


//*********************************************************************************************

/*! @brief Compute energytransfer(ET) nlin \f$ N_i \leftarrow D_j \mathcal{F}(Wr_j * G_i) \f$
 *			G is the giver field (G can be V or W).
 *
 *  @param	W IncVF that acts as a helper field.
 *
 *	Same as void IncVF::EnergyTr_Compute_nlin() except that here W plays the role of helper
 *		in place of V.
 *
 *  @return \f$ N_i \leftarrow D_j \mathcal{F}(Wr_j * G_i) \f$.
 */
void IncVF::EnergyTr_Compute_nlin(IncVF& W) 
{

#ifndef TRANSPOSE
	*W.V1r = *W.V1;  
	*W.V2r = *W.V2; 
	*W.V3r = *W.V3;    
	
	W.RV_Inverse_transform(*VF_temp_r);
	
	ET_Inverse_transform_G(*VF_temp_r);

	Compute_RSprod_diagET(W);			
 
	NLIN_diag_Forward_transform_derivative(*VF_temp_r);			
	  
	Compute_RSprod_offdiagET(W);		
	
	RV_Forward_transform_RSprod(*VF_temp_r);				// Vr[i]=FFT(Vr[i]*Vr[j])
	ET_Forward_transform_RSprod_G(*VF_temp_r);
	
	Derivative_RSprod_VV_ET();								
	// Compute derivatives and construct nlin
	
//
//
	
#else
	W.RV_Inverse_transform_transpose_order(*W.V1, *W.V2, *W.V3, *W.VF_temp);  
	
	ET_Inverse_transform_G_transpose_order();
					                      
	Compute_RSprod_diag_ft_derivativeET(W);			
 
	Compute_RSprod_offdiagET(W);					
	
	Forward_transform_derivative_RSprod_offdiagET();	
	
#endif



}


//*********************************************************************************************

/*! @brief Compute energytransfer(ET) nlin \f$ N_i \leftarrow D_j \mathcal{F}(Vr_j * T^<) \f$
 *			\f$T^< \f$ is the giver field.  \f$T^< \f$ is kept in G1 array.
 *
 *	@sa		ET_compute_nlin.cc
 *
 *  @return \f$ N_i \leftarrow D_j \mathcal{F}(Vr_j * T^<) \f$.
 */
void IncVF::EnergyTr_Compute_nlin(IncSF& T) 
{
  
#ifndef TRANSPOSE

	*V1r = *V1;  
	*V2r = *V2;  
	*V3r = *V3;          
	
	RVF::RV_Inverse_transform(*VF_temp_r);
	
	Inverse_transform_array(basis_type, N, *G1, *VF_temp_r, 1);		// *G = *T.F<; parity =1 (sin)
	
	Array_real_mult(N, *V1r, *G1, *nlin1);			
	Array_real_mult(N, *V2r, *G1, *nlin2);
	Array_real_mult(N, *V3r, *G1, *nlin3);		
	
	Forward_transform_array(basis_type, N, *nlin1, *VF_temp_r, 0);	
	// parity = 0; sin*sin = cos	
	
	Forward_transform_array(basis_type, N, *nlin2, *VF_temp_r, 1);		
	// parity = 1; sin*cos = sin
	
	Forward_transform_array(basis_type, N, *nlin3, *VF_temp_r, 1);		
	// parity = 1; sin*cos = sin
	
	Xderiv(basis_type, N, *nlin1, *nlin1, kfactor, 0);
	Yderiv(basis_type, N, *nlin2, *nlin2, kfactor, 0);
	Zderiv(basis_type, N, *nlin3, *nlin3, kfactor, 0);

	*nlin1 = *nlin1 + *nlin2;	
	*nlin1 = *nlin1 + *nlin3;	

//
//

#else

	RVF::RV_Inverse_transform_transpose_order(*V1, *V2, *V3, *VF_temp);    
	
	// *G = *T.F<; parity =1 (sin)
	*VF_temp_r = *G1;
	Inverse_transform_array_transpose_order(basis_type, N, *VF_temp_r, *G1, 1);		
	//  Now G1 is in real space (in transpose order)					

	Array_real_mult(N, *V1r, *G1, *VF_temp_r);
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin1, 0);
	Xderiv(basis_type, N, *nlin1, *nlin1, kfactor, 0);
									
	Array_real_mult(N, *V2r, *G1, *VF_temp_r);			
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin2, 1);
	Yderiv(basis_type, N, *nlin2, *nlin2, kfactor, 0);
	
	Array_real_mult(N, *V3r, *G1, *VF_temp_r);			
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin3, 1);
	Zderiv(basis_type, N, *nlin3, *nlin3, kfactor, 0);
	
	*nlin1 = *nlin1 + *nlin2;	
	*nlin1 = *nlin1 + *nlin3;
	
#endif	

}



//*********************************************************************************************

/*! @brief Compute energytransfer(ET) nlin \f$ N_i \leftarrow D_j \mathcal{F}(\omega_j * G_i) \f$
 *			G is the giver field (G can be V or W).
 *
 *	@sa ET_compute_nlin.cc
 *
 *  @return \f$ N_i \leftarrow D_j \mathcal{F}(\omega_j * V_i) \f$.
 */ 
void IncVF::EnergyTr_Compute_nlin_vorticity_helper()
{  
	
	TinyVector<complx,3> vorticity;
	
	for (int i1=0;  i1<N[1]; i1++)
		for (int i2=0;  i2<N[2]; i2++) 
			for (int i3=0;  i3<=N[3]/2; i3++) 
			{
				if (N[2] > 1)
				{	
					Compute_Modal_vorticity(basis_type, i1, i2, i3, N, 
											*V1, *V2, *V3, kfactor, vorticity);
					
					(*V1r)(i1, i2, i3) = vorticity(0);
					(*V2r)(i1, i2, i3) = vorticity(1);
					(*V3r)(i1, i2, i3) = vorticity(2);
				}	
				
				else if (N[2] == 1)
				{
					cout << " ERROR: In EnergyTr_Compute_nlin_vorticity_helper() " << endl;
					cout << "N[2] must be greater than 1 for hm calculations " << endl;
					exit(0);
				}
			}
	
	
	RV_Inverse_transform(*VF_temp_r);		// Vr would contain omega(r)
	ET_Inverse_transform_G(*VF_temp_r);
	
	IncVF::Compute_RSprod_diagET();					
	
	IncVF::NLIN_diag_Forward_transform_derivative(*VF_temp_r);                     
	
	IncVF::Compute_RSprod_offdiagET();	
	
	RVF::RV_Forward_transform_RSprod(*VF_temp_r);					// V_i = T[ Hj G_i] 
	
	ET_Forward_transform_RSprod_G(*VF_temp_r);
	
	IncVF::Derivative_RSprod_VV_ET();	 	
}

	//*********************************************************************************************

/*! @brief Compute energytransfer(ET) nlin \f$ N_i Vr_j \times G_i) \f$
 *			G is the giver field  W.
 *
 *	@sa ET_compute_nlin.cc
 *
 *  @return \f$ N_i \leftarrow D_j \mathcal{F}(Vr_j * V_i) \f$.
 */ 
void IncVF::EnergyTr_Compute_nlin_UcrossB()
{  
	
	*V1r = *V1;  
	*V2r = *V2; 
	*V3r = *V3;   	
	
	RV_Inverse_transform(*VF_temp_r);		
	
	ET_Inverse_transform_G(*VF_temp_r);
	
	if (N[2] > 1)
	{	
		*nlin1 = (*V2r)*(*G3) - (*V3r)*(*G2);
		*nlin2 = (*V3r)*(*G1) - (*V1r)*(*G3);
		*nlin3 = (*V1r)*(*G2) - (*V2r)*(*G1);
		
		Forward_transform_array(basis_type, N, *nlin1, *VF_temp_r, 0);
		Forward_transform_array(basis_type, N, *nlin2, *VF_temp_r, 0);
		Forward_transform_array(basis_type, N, *nlin3, *VF_temp_r, 0);
	}
	
	else if (N[2] == 1)
	{	
		*nlin1 = (*V2r)*(*G3) - (*V3r)*(*G2);
		*nlin3 = (*V1r)*(*G2) - (*V2r)*(*G1);
		
		Forward_transform_array(basis_type, N, *nlin1, *VF_temp_r, 0);
		Forward_transform_array(basis_type, N, *nlin3, *VF_temp_r, 0);
	}
	
}


//*****************************  End of ET_compute_nlin.cc  ***********************************





