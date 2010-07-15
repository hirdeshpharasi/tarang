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
 * @brief  Compute the diagonal and nondiagonal terms of the real space products 
 *			(TRANSPOSE ORDER).
 *
 * @note   Assume N1 = N2 (VF_Temp = VF_temp_r)
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug   TO BE TESTED
 */

#include "../IncVF.h"
#include "../IncSF.h"


//*********************************************************************************************

/*! @brief Compute real space diagonal product \f$ (Vr_i)^2 \f$ and 
 *			\f$ N_i \leftarrow D_i \mathcal{F}(Vr_i)^2 \f$ [no i sum].
 *
 *  For fluid (a) \f$ N_i = (Vr_i)^2 \f$ (real space). 
 *            (b) Forward transform $N_i$ in transpose order
 *            (c) Compute \f$ N_i \leftarrow D_i \mathcal{F}(Vr_i)^2  \f$.
 *
 * @return \f$ N_i \leftarrow D_i \mathcal{F}(Vr_i * Wr_i) \f$.
 */
void IncVF::Compute_RSprod_diag_ft_derivative()
{
#ifdef TRANSPOSE
	Array_real_mult(N, *V1r, *V1r, *VF_temp_r);
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin1, 0);
	Xderiv(basis_type, N, *nlin1, *nlin1, kfactor, 0);
	
	Array_real_mult(N, *V2r, *V2r, *VF_temp_r);
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin2, 0);
	Yderiv(basis_type, N, *nlin2, *nlin2, kfactor, 0);  // parity unimportant here.
	
	Array_real_mult(N, *V3r, *V3r, *VF_temp_r); 
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin3, 0);
	Zderiv(basis_type, N, *nlin3, *nlin3, kfactor, 0);	 // parity unimportant here.
#endif		
}


//*********************************************************************************************

/*! @brief Compute real space diagonal product \f$ (Vr_i)^2 \f$ and 
 *			\f$ N_i \leftarrow D_i \mathcal{F}(Vr_i * Wr_i) \f$ [no i sum].
 *
 *  For fluid (a) \f$ N_i = (Vr_i Wr_i) \f$ (real space). 
 *            (b) Forward transform $N_i$ in transpose order
 *            (c) Compute \f$ N_i \leftarrow D_i \mathcal{F}(Vr_i * Wr_i) \f$.
 *
 * @return \f$ N_i \leftarrow D_i \mathcal{F}(Vr_i * Wr_i) \f$.
 */
void IncVF::Compute_RSprod_diag_ft_derivative(IncVF& W)
{
#ifdef TRANSPOSE
	Array_real_mult(N, *V1r, *W.V1r, *VF_temp_r);
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin1, 0);
	Xderiv(basis_type, N, *nlin1, *nlin1, kfactor, 0);  // parity unimportant here.
	
	Array_real_mult(N, *V2r, *W.V2r, *VF_temp_r);
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin2, 0);
	Yderiv(basis_type, N, *nlin2, *nlin2, kfactor, 0);  // parity unimportant here.
	
	Array_real_mult(N, *V3r, *W.V3r, *VF_temp_r); 
	Forward_transform_array_transpose_order(basis_type, N, *VF_temp_r, *nlin3, 0);
	Zderiv(basis_type, N, *nlin3, *nlin3, kfactor, 0);	// parity unimportant here.
#endif	 	
}

//*********************************************************************************************
/*! @brief Compute real space off-diagonal product \f$ Vr_i*Vr_j \f$.
 *
 *  Before entering the function <BR>  \f$ Vr_i = (Vr_i * Vr_j) \f$. <BR>
 *
 * Steps:  (a) Compute forward transform in transpose order.
 *         (b) Compute derivatives.
 *
 *  @return \f$ N_i \leftarrow D_j \mathcal{F}(Vr_i * Vr_j) \f$.
 */
void IncVF::Forward_transform_derivative_RSprod_offdiag()
{
#ifdef TRANSPOSE
	Forward_transform_array_transpose_order(basis_type, N, *V1r, *VF_temp, 1);	// SFT
	Yderiv_RSprod_VV(*VF_temp, *VF_temp2);   *nlin1 = *nlin1 + *VF_temp2; 
	Xderiv_RSprod_VV(*VF_temp, *VF_temp2);   *nlin2 = *nlin2 + *VF_temp2; 
	
	Forward_transform_array_transpose_order(basis_type, N, *V2r, *VF_temp, 0);	// CFT
	Zderiv_RSprod_VV(*VF_temp, *VF_temp2);   *nlin2 = *nlin2 + *VF_temp2;
	Yderiv_RSprod_VV(*VF_temp, *VF_temp2);   *nlin3 = *nlin3 + *VF_temp2; 
	
	Forward_transform_array_transpose_order(basis_type, N, *V3r, *VF_temp, 1);	// SFT
	Zderiv_RSprod_VV(*VF_temp, *VF_temp2);   *nlin1 = *nlin1 + *VF_temp2;
	Xderiv_RSprod_VV(*VF_temp, *VF_temp2);   *nlin3 = *nlin3 + *VF_temp2; 
#endif
}



//*********************************************************************************************
/*! @brief Compute real space off-diagonal product \f$ Vr_i*Vr_j \f$.
 *
 *  Before entering the function <BR>  \f$ Vr_i = (Vr_i * Wr_j) \f$. <BR>
 *
 * Steps:  (a) Compute forward transform in transpose order.
 *         (b) Compute derivatives.
 *
 *  @return \f$ N_i \leftarrow D_j \mathcal{F}(Vr_i * Wr_j) \f$.
 */
void IncVF::Forward_transform_derivative_RSprod_offdiag(IncVF& W)
{
#ifdef TRANSPOSE
	Forward_transform_array_transpose_order(basis_type, N, *V1r, *VF_temp, 1);		// SFT
	Yderiv_RSprod_VV(*VF_temp, *VF_temp2);   *nlin1 = *nlin1 + *VF_temp2; 
	Xderiv_RSprod_VV(*VF_temp, *VF_temp2);   *W.nlin2 = *W.nlin2 + *VF_tem2p; 
	
	Forward_transform_array_transpose_order(basis_type, N, *V2r, *VF_temp, 0);		// CFT
	Zderiv_RSprod_VV(*VF_temp, *VF_temp2);   *nlin2 = *nlin2 + *VF_temp2;
	Yderiv_RSprod_VV(*VF_temp, *VF_temp2);   *W.nlin3 = *W.nlin3 + *VF_temp2; 
	
	Forward_transform_array_transpose_order(basis_type, N, *V3r, *VF_temp, 1);		// SFT
	Zderiv_RSprod_VV(*VF_temp, *VF_temp2);   *W.nlin1 = *W.nlin1 + *VF_temp2;
	Xderiv_RSprod_VV(*VF_temp, *VF_temp2);   *nlin3 = *nlin3 + *VF_temp2; 
		
	Forward_transform_array_transpose_order(basis_type, N, *W.V1r, *VF_temp, 1);	// SFT
	Yderiv_RSprod_VV(*VF_temp, *VF_temp2);   *W.nlin1 = *W.nlin1 + *VF_temp2; 
	Xderiv_RSprod_VV(*VF_temp, *VF_temp2);   *nlin2 = *nlin2 + *VF_temp2; 
	
	Forward_transform_array_transpose_order(basis_type, N, *W.V2r, *VF_temp, 0);	// CFT
	Zderiv_RSprod_VV(*VF_temp, *VF_temp2);   *W.nlin2 = *W.nlin2 + *VF_temp2;
	Yderiv_RSprod_VV(*VF_temp, *VF_temp2);   *nlin3 = *nlin3 + *VF_temp2; 
	
	Forward_transform_array_transpose_order(basis_type, N, *W.V3r, *VF_temp, 1);	// SFT
	Zderiv_RSprod_VV(*VF_temp, *VF_temp2);   *nlin1 = *nlin1 + *VF_temp2;
	Xderiv_RSprod_VV(*VF_temp, *VF_temp2);   *W.nlin3 = *W.nlin3 + *VF_temp2; 
#endif
}

//****************************  End of RSprod_transpose.cc ************************************


