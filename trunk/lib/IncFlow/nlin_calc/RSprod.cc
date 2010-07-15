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


/*! \file  RSprod.cc
 * 
 * @brief  Compute the diagonal and nondiagonal terms of the real space products 
 *			(NORMAL ORDER).
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug	No known bugs
 */


#include "../IncVF.h"
#include "../IncSF.h"


//*********************************************************************************************


/*! @brief Compute real space diagonal product \f$ (Vr_i)^2 \f$.
 *
 *  @return  \f$ N_i = (Vr_i)^2 \f$. 
 */
void IncVF::Compute_RSprod_diag()
{
	Array_real_mult(N, *V1r, *V1r, *nlin1);  
	Array_real_mult(N, *V2r, *V2r, *nlin2);
	Array_real_mult(N, *V3r, *V3r, *nlin3);   
}



//*********************************************************************************************


/*! @brief Compute real space off-diagonal product \f$ Vr_i*Vr_j \f$.
 *
 *  @return 3D: \f$ Vr_1 = (Vr_1 * Vr_2) \f$; <BR> 
 *				\f$ Vr_2 = (Vr_2 * Vr_3) \f$; <BR>	 
 *				\f$ Vr_3 = (Vr_3 * Vr_1) \f$.
 */
void IncVF::Compute_RSprod_offdiag()
{
	*VF_temp = *V1r;  

	Array_real_mult(N, *V1r, *V2r, *V1r);		// V1r= V1r*V2r
	Array_real_mult(N, *V2r, *V3r, *V2r);		// V2r= V2r*V3r
	Array_real_mult(N, *V3r, *VF_temp, *V3r);	// V3r= V3r*V1r 
}



//*********************************************************************************************

/*! @brief Compute \f$ N_i \leftarrow D_j \mathcal{F}(Vr_i * Vr_j) \f$.
 *
 *	Before entering the function <BR>
 *		\f$ N_i = D_i \mathcal{F}(Vr_i^2) \f$ <BR>
 *
 *		3D: \f$ Vr_1 = \mathcal{F}(Vr_1 * Vr_2) \f$; <BR>
 *			\f$ Vr_2 = \mathcal{F}(Vr_2 * Vr_3) \f$; <BR>
 *			\f$	Vr_3 = \mathcal{F}(Vr_3 * Vr_1) \f$.
 *
 *  @return \f$ N_i \leftarrow D_j \mathcal{F}(Vr_i * Vr_j) \f$.
 */
void IncVF::Derivative_RSprod_VV()  
{
  
	Yderiv_RSprod_VV(*V1r, *VF_temp);   *nlin1 = *nlin1 + *VF_temp; 
	Zderiv_RSprod_VV(*V3r, *VF_temp);   *nlin1 = *nlin1 + *VF_temp;       
 
	Xderiv_RSprod_VV(*V1r, *VF_temp);   *nlin2 = *nlin2 + *VF_temp; 
	Zderiv_RSprod_VV(*V2r, *VF_temp);   *nlin2 = *nlin2 + *VF_temp;      
	
	Xderiv_RSprod_VV(*V3r, *VF_temp);   *nlin3 = *nlin3 + *VF_temp;   
	Yderiv_RSprod_VV(*V2r, *VF_temp);   *nlin3 = *nlin3 + *VF_temp;           
}



//*********************************************************************************************

/*! @brief Compute real space diagonal product \f$ (Vr_i*Wr_i) \f$.
 *
 *	@param	W  IncVF, e.g, for magnetic field B for MHD
 *
 *  @return  \f$ N_i = (Vr_i * Wr_i) \f$. 
 */
void IncVF::Compute_RSprod_diag(IncVF& W)           
{
	Array_real_mult(N, *V1r, *W.V1r, *nlin1);  
	Array_real_mult(N, *V2r, *W.V2r, *nlin2);
	Array_real_mult(N, *V3r, *W.V3r, *nlin3);   
}



//*********************************************************************************************

/*! @brief Compute real space off-diagonal product \f$ Vr_i*Vr_j \f$.
 *
 *  @return 3D: \f$ Vr_1 = (Vr_1 * Wr_2) \f$; <BR> 
 *				\f$ Vr_2 = (Vr_2 * Wr_3) \f$; <BR>	 
 *				\f$ Vr_3 = (Vr_3 * Wr_1) \f$. <BR>
 *
 *				\f$ Wr_1 = (Wr_1 * Vr_2) \f$; <BR> 
 *				\f$ Wr_2 = (Wr_2 * Vr_3) \f$; <BR>	 
 *				\f$ Wr_3 = (Wr_3 * Vr_1) \f$. <BR>
 */
void IncVF::Compute_RSprod_offdiag(IncVF& W)
{
	*VF_temp = *V1r;  
	*W.VF_temp = *W.V1r;								// store V1r and W1r in VF_temp 
  
	Array_real_mult(N, *V1r,   *W.V2r,  *V1r);			// V1r = V1r*W2r
	Array_real_mult(N, *W.V1r, *V2r,    *W.V1r);		// W1r = W1r*V2r
	
	Array_real_mult(N, *V2r,   *W.V3r,  *V2r);			// V2r = V2r*W3r
	Array_real_mult(N, *W.V2r, *V3r,     *W.V2r);		// W2r = W2r*V3
	
	Array_real_mult(N, *V3r,   *W.VF_temp, *V3r);		// V3r = V3r*W1r
	Array_real_mult(N, *W.V3r, *VF_temp,   *W.V3r);		// W3r = W3r*V1r		
}


//*********************************************************************************************

/*! @brief Compute \f$ N_i \leftarrow D_j \mathcal{F}(Vr_i * Wr_j) \f$ and
 *				   \f$ N^W_i \leftarrow D_j \mathcal{F}(Wr_i * Vr_j) \f$.
 *
 *	Before entering the function <BR>
 *			\f$ N_i   = D_i \mathcal{F}(Vr_i * Wr_i) \f$ [no i sum] and <BR> 
 *			\f$ N^W_i = D_i \mathcal{F}(Vr_i * Wr_i) \f$ <BR>
 *
 *		3D: \f$ Vr_1 = \mathcal{F}(Vr_1 * Wr_2) \f$; <BR>
 *			\f$ Vr_2 = \mathcal{F}(Vr_2 * Wr_3) \f$; <BR> 
 *			\f$	Vr_3 = \mathcal{F}(Vr_3 * Wr_1) \f$, and <BR>
 *
 *			\f$ Wr_1 = \mathcal{F}(Wr_1 * Vr_2) \f$; <BR>
 *			\f$ Wr_2 = \mathcal{F}(Wr_2 * Vr_3) \f$; <BR>
 *			\f$	Wr_3 = \mathcal{F}(Wr_3 * Vr_1) \f$
 *
 *  @return \f$ N_i \leftarrow D_j \mathcal{F}(Vr_i * Wr_j) \f$ and
 *				   \f$ N^W_i \leftarrow D_j \mathcal{F}(Wr_i * Vr_j) \f$.
 */
void IncVF::Derivative_RSprod_VV(IncVF& W)  
{
	Yderiv_RSprod_VV(*V1r, *VF_temp);		*nlin1 = *nlin1 + *VF_temp; 
	Zderiv_RSprod_VV(*W.V3r, *VF_temp);		*nlin1 = *nlin1 + *VF_temp; 
	
	Xderiv_RSprod_VV(*W.V1r, *VF_temp);		*nlin2 = *nlin2 + *VF_temp; 
	Zderiv_RSprod_VV(*V2r, *VF_temp);		*nlin2 = *nlin2 + *VF_temp;   
  
  	Xderiv_RSprod_VV(*V3r, *VF_temp);		*nlin3 = *nlin3 + *VF_temp;   
	Yderiv_RSprod_VV(*W.V2r, *VF_temp);		*nlin3 = *nlin3 + *VF_temp; 
	
	Yderiv_RSprod_VV(*W.V1r, *VF_temp);		*W.nlin1 = *W.nlin1 + *VF_temp; 
	Zderiv_RSprod_VV(*V3r, *VF_temp);		*W.nlin1 = *W.nlin1 + *VF_temp;  
	
	Xderiv_RSprod_VV(*V1r, *VF_temp);		*W.nlin2 = *W.nlin2 + *VF_temp; 
	Zderiv_RSprod_VV(*W.V2r, *VF_temp);		*W.nlin2 = *W.nlin2 + *VF_temp;      
	
	Xderiv_RSprod_VV(*W.V3r, *VF_temp);		*W.nlin3 = *W.nlin3 + *VF_temp;   
	Yderiv_RSprod_VV(*V2r, *VF_temp);		*W.nlin3 = *W.nlin3 + *VF_temp;    

}


//****************************  End of RSprod.cc **********************************************



