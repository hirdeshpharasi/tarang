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

/*! \file  ET_RSprod_transpose.cc
 * 
 * @brief  Compute the diagonal and nondiagonal terms of the real space products
 *			(Transpose order)
 *
 * @note  Assume N1 = N2 for transpose order.  The dimension of G_i(real space) is assume
 *			to be the same as G_i(Fourie space).
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

/*! @brief Compute real space diagonal product  \f$ (Vr_i * G_i) \f$ and 
 *			\f$ N_i \leftarrow D_i \mathcal{F}(Vr_i * G_i) \f$ [no i sum].
 *
 *  @note   Here Vi and Gi while entering are in transpose order (real arrays).
 *
 *  Steps	  (a) \f$ N_i = (Vr_i * G_i) \f$ (real space). 
 *            (b) Forward transform $N_i$ in transpose order
 *            (c) Compute \f$ N_i \leftarrow D_i \mathcal{F}(Vr_i * G_i)  \f$.
 *
 * @return \f$ N_i \leftarrow D_i \mathcal{F}(Vr_i * G_i) \f$.
 */
void IncVF::Compute_RSprod_diag_ft_derivativeET()
{
#ifdef TRANSPOSE 
	Array_real_mult(N, *V1r, *G1, *VF_temp_r);    
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin1, 0);
	Xderiv(basis_type, N, *nlin1, *nlin1, kfactor, 0);
	
	Array_real_mult(N, *V2r, *G2, *VF_temp_r);
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin2, 0);
	Yderiv(basis_type, N, *nlin2, *nlin2, kfactor, 0);
	
	Array_real_mult(N, *V3r, *G3, *VF_temp_r); 
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin3, 0);
	Zderiv(basis_type, N, *nlin3, *nlin3, kfactor, 0);	 	
#endif	
}


//*********************************************************************************************
/*! @brief Compute real space diagonal product  \f$ (Wr_i * G_i) \f$ and 
 *			\f$ N_i \leftarrow D_i \mathcal{F}(Wr_i * G_i) \f$ [no i sum].
 *
 *  @note   Here Vi and Gi while entering are in transpose order (real arrays).
 *
 *  Steps	  (a) \f$ N_i = (Wr_i * G_i) \f$ (real space). 
 *            (b) Forward transform $N_i$ in transpose order
 *            (c) Compute \f$ N_i \leftarrow D_i \mathcal{F}(Wr_i * G_i)  \f$.
 *
 * @return \f$ N_i \leftarrow D_i \mathcal{F}(Wr_i * G_i) \f$.
 */
void IncVF::Compute_RSprod_diag_ft_derivativeET(IncVF& W)
{
#ifdef TRANSPOSE 
	Array_real_mult(N, *W.V1r, *G1, *VF_temp_r);
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin1, 0);
	Xderiv(basis_type, N, *nlin1, *nlin1, kfactor, 0);
	
	Array_real_mult(N, *W.V2r, *G2, *VF_temp_r);
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin2, 0);
	Yderiv(basis_type, N, *nlin2, *nlin2, kfactor, 0);
	
	Array_real_mult(N, *W.V3r, *G3, *VF_temp_r); 
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin3, 0);
	Zderiv(basis_type, N, *nlin3, *nlin3, kfactor, 0);	 	
#endif	
}



//*********************************************************************************************
/*! @brief Compute real space diagonal product  \f$ (Vr_i * G_j) \f$ and 
 *			\f$ N_i \leftarrow D_i \mathcal{F}(Vr_i * G_i) \f$ [no i sum].
 *
 *  @note   Here Vi and Gi while entering are in transpose order (real arrays).
 *
 *  Steps	  (a) Compute real space diagonal product  \f$ (Vr_i * G_j) \f$.
 *            (b) Forward transform the products in transpose order
 *            (c) Compute \f$ N_i \leftarrow D_i \mathcal{F}\f$ (Vr_i * G_j)  \f$.
 *
 * @return \f$ N_i \leftarrow D_i \mathcal{F}\f$ (Vr_i * G_j) \f$.
 */
void IncVF::Forward_transform_derivative_RSprod_offdiagET()
{
#ifdef TRANSPOSE
	
	Forward_transform_array_transpose_order(basis_type, N, *V1r, *VF_temp, 1);	
	Xderiv_RSprod_VV(*VF_temp, *VF_temp);   
	*nlin2 = *nlin2 + *VF_temp; 
	
	Forward_transform_array_transpose_order(basis_type, N, *V2r, *VF_temp, 0);	
	Yderiv_RSprod_VV(*VF_temp, *VF_temp);  
	*nlin3 = *nlin3 + *VF_temp; 
	
	Forward_transform_array_transpose_order(basis_type, N, *V3r, *VF_temp, 1);
	Zderiv_RSprod_VV(*VF_temp, *VF_temp);   
	*nlin1 = *nlin1 + *VF_temp;
	
	Forward_transform_array_transpose_order(basis_type, N, *G1, *VF_temp, 1);	
	Yderiv_RSprod_VV(*VF_temp, *VF_temp);   
	*nlin1 = *nlin1 + *VF_temp; 
	
	Forward_transform_array_transpose_order(basis_type, N, *G2, *VF_temp, 0);	
	Zderiv_RSprod_VV(*VF_temp, *VF_temp);   
	*nlin2 = *nlin2 + *VF_temp; 
	
	Forward_transform_array_transpose_order(basis_type, N, *G3, *VF_temp, 1);
	Xderiv_RSprod_VV(*VF_temp, *VF_temp);   
	*nlin3 = *nlin3 + *VF_temp;
#endif	
}

//********************************  End of ET_RSprod_transform.cc *****************************








