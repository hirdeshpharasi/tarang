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

/*! \file  ET_RSprod.cc
 * 
 * @brief  Compute the diagonal and nondiagonal terms of the real space products.
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

/*! @brief Compute real space diagonal product \f$ (Vr_i * G_i) \f$.
 *
 *  @return  \f$ N_i \leftarrow (Vr_i * G_i) \f$. 
 */
 void IncVF::Compute_RSprod_diagET()           
{
	Array_real_mult(N, *V1r, *G1, *nlin1);  
	Array_real_mult(N, *V2r, *G2, *nlin2);
	Array_real_mult(N, *V3r, *G3, *nlin3);   
}


//*********************************************************************************************

/*! @brief Compute real space off-diagonal product \f$ Vr_i*G_j \f$.
 *
 *  @return 2D: \f$ Vr_1 = (Vr_1 * G_2) \f$. <BR>
 *				\f$ G_1 = (G_1 * Vr_2) \f$. 
 *
 *  @return 3D: \f$ Vr_1 = (Vr_1 * G_2) \f$; <BR> 
 *				\f$ Vr_2 = (Vr_2 * G_3) \f$; <BR>	 
 *				\f$ Vr_3 = (Vr_3 * G_1) \f$. <BR>
 *
 *				\f$ G_1 = (G_1 * Vr_2) \f$; <BR> 
 *				\f$ G_2 = (G_2 * Vr_3) \f$; <BR>	 
 *				\f$ G_3 = (G_3 * Vr_1) \f$.
 */
void IncVF::Compute_RSprod_offdiagET()
{
	*VF_temp = *V1r;    
	*VF_temp2 = *G1;				// store V1r and W1r in temp 
  
	Array_real_mult(N, *V1r, *G2,	*V1r);			// V1r = V1r*G2
	Array_real_mult(N, *G1,  *V2r,	*G1);			// G1 = G1*V2r
	
	Array_real_mult(N, *V2r, *G3,	*V2r);			// V2r = V2r*G3
	Array_real_mult(N, *G2,	 *V3r,  *G2);			// G2 =  G2*V3r
	
	Array_real_mult(N, *V3r, *VF_temp2, *V3r);		// V3r = V3r*G1
	Array_real_mult(N, *G3,  *VF_temp, *G3);		// G3 = G3*V1r
}



//*********************************************************************************************

/*! @brief Compute real space diagonal product \f$ (Wr_i * G_i) \f$.
 *
 *	@param	W  IncVF, e.g., magnetic field for MHD.
 *
 *  @return  \f$ N_i \leftarrow (Wr_i * G_i) \f$. 
 */
void IncVF::Compute_RSprod_diagET(IncVF& W)           
{
	Array_real_mult(N, *W.V1r, *G1, *nlin1);  
	Array_real_mult(N, *W.V2r, *G2, *nlin2);
	Array_real_mult(N, *W.V3r, *G3, *nlin3);   
}


//*********************************************************************************************

/*! @brief Compute real space off-diagonal product \f$ Wr_i*Gr_j \f$.
 *
 *	@param	W  IncVF, e.g., magnetic field for MHD.
 *
 *  @return 2D: \f$ Vr_1 = (Wr_1 * G_2) \f$. <BR>
 *				\f$ G_1  = (G_1 * Wr_2) \f$. 
 *
 *  @return 3D: \f$ Vr_1 = (Wr_1 * G_2) \f$; <BR> 
 *				\f$ Vr_2 = (Wr_2 * G_3) \f$; <BR>	 
 *				\f$ Vr_3 = (Wr_3 * G_1) \f$. <BR>
 *
 *				\f$ G_1 = (G_1 * Wr_2) \f$; <BR> 
 *				\f$ G_2 = (G_2 * Wr_3) \f$; <BR>	 
 *				\f$ G_3 = (G_3 * Wr_1) \f$.
 */
void IncVF::Compute_RSprod_offdiagET(IncVF& W)
{
  
	Array_real_mult(N, *W.V1r, *G2, *V1r);			// V1r = W1r*G2
	Array_real_mult(N, *W.V2r, *G3, *V2r);			// V2r = W2r*G3
	Array_real_mult(N, *W.V3r, *G1, *V3r);			// V3r = W3r*G1
	
	Array_real_mult(N, *G1, *W.V2r,	 *G1);			// G1 = G1*W2r	
	Array_real_mult(N, *G2, *W.V3r,	 *G2);			// G2 = G2*W3r
	Array_real_mult(N, *G3, *W.V1r,  *G3);			// G3 = G3*W1r
}



//*********************************************************************************************

/*! @brief Compute \f$ N_i \leftarrow D_j \mathcal{F}(G_i * Vr_j) \f$.
 *			Works for both \f$ (U.\nabla) G \f$ and \f$ (W.\nabla) G \f$.
 *
 *	Before entering the function <BR>
 *			\f$ N_i = D_i \mathcal{F}(Vr_i * G_i) \f$.
 *
 *		2D: \f$ Vr_1 = \mathcal{F}(Vr_1 * G_2) \f$; <BR>
 *			\f$ G_1 = \mathcal{F}(G_1 * Vr_2) \f$. <BR>
 *
 *		3D: \f$ Vr_1 = \mathcal{F}(Vr_1 * G_2) \f$; <BR>
 *			\f$ Vr_2 = \mathcal{F}(Vr_2 * G_3) \f$; <BR> 
 *			\f$	Vr_3 = \mathcal{F}(Vr_3 * G_1) \f$, and <BR>
 *
 *			\f$ G_1 = \mathcal{F}(G_1 * Vr_2) \f$; <BR>
 *			\f$ G_2 = \mathcal{F}(G_2 * Vr_3) \f$; <BR>
 *			\f$	G_3 = \mathcal{F}(G_3 * Vr_1) \f$
 *
 *  @return \f$ N_i \leftarrow D_j \mathcal{F}(Vr_i * Wr_j) \f$ and
 *				   \f$ N^W_i \leftarrow D_j \mathcal{F}(Wr_i * Vr_j) \f$.
 */
void IncVF::Derivative_RSprod_VV_ET()  
{

	Yderiv_RSprod_VV(*G1, *VF_temp);		*nlin1 = *nlin1 + *VF_temp; 
	Zderiv_RSprod_VV(*V3r, *VF_temp);		*nlin1 = *nlin1 + *VF_temp;  
	
	Xderiv_RSprod_VV(*V1r, *VF_temp);		*nlin2 = *nlin2 + *VF_temp; 
	Zderiv_RSprod_VV(*G2, *VF_temp);		*nlin2 = *nlin2 + *VF_temp;       
	
	Xderiv_RSprod_VV(*G3, *VF_temp);		*nlin3 = *nlin3 + *VF_temp;   
	Yderiv_RSprod_VV(*V2r, *VF_temp);		*nlin3 = *nlin3 + *VF_temp;
}


//****************************  End of ET_RSprod.cc *******************************************






