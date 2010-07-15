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

/*! \file  compute_derivative.cc
 * 
 * @brief  Compute the diagonal and nondiagonal terms of the real space products.
 *
 *	Before entering the functions deriv_RSprod_VV <BR>
 *		\f$ N_i = \mathcal{F}(Vr_i^2) \f$ <BR>
 *
 *			\f$ Vr_1 = \mathcal{F}(Vr_1 * Vr_2) \f$; <BR>
 *			\f$ Vr_2 = \mathcal{F}(Vr_2 * Vr_3) \f$; <BR>
 *			\f$	Vr_3 = \mathcal{F}(Vr_3 * Vr_1) \f$.
 *
 *	Before entering the functions deriv_RSprod_VT <BR>
 *		\f$ N^T_i = \mathcal{F}(Vr_i * Fr) \f$ <BR>
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec 2008
 *
 * @bug		No known bugs
 */

#include "../IncVF.h"
#include "../IncSF.h"


//*********************************************************************************************

// For Real space product Vi x Vj 

/*! @brief Compute \f$ B \leftarrow D_x(A) \f$.
 *
 *	@note Before entering this function, \f$ A = Vr_x Vr_j \f$ with j=2 or 3. <BR>
 *			Thus A has odd parity (sin*cos = sin).
 *
 *  @return \f$ B \leftarrow D_x(A) \f$.
 */
void IncVF::Xderiv_RSprod_VV(Array<complx,3> A, Array<complx,3> B)
{
	Xderiv(basis_type, N, A, B, kfactor, 1);		// Parity = 1; sin*cos = sin
}	


//
//

/*! @brief Compute \f$ B \leftarrow D_y(A) \f$.
 *
 *	@note The basis functions are Fourier along y direction. So parity does not matter. 
 *
 *  @return \f$ B \leftarrow D_y(A) \f$.
 */
void IncVF::Yderiv_RSprod_VV(Array<complx,3> A, Array<complx,3> B)
{
	Yderiv(basis_type, N, A, B, kfactor, 0);
}


// Zderive
	
/*! @brief Compute \f$ B \leftarrow D_z(A) \f$.
 *
 *	@note The basis functions are Fourier along z direction. So parity does not matter. 
 *
 *  @return \f$ B \leftarrow D_z(A) \f$.
 */
void IncVF::Zderiv_RSprod_VV(Array<complx,3> A, Array<complx,3> B)
{
	Zderiv(basis_type, N, A, B, kfactor, 0);
}


//*********************************************************************************************


/*! @brief Compute \f$ B \leftarrow D_x(A) \f$.
 *
 *	@note Before entering this function, \f$ A = T Vr_x \f$ <BR>
 *			Thus A has even parity (sin*sin = cos).
 *
 *  @return \f$ B \leftarrow D_x(A) \f$.
 */
void IncVF::Xderiv_RSprod_VT(Array<complx,3> A, Array<complx,3> B)
{
	Xderiv(basis_type, N, A, B, kfactor, 0);		// Parity = 0; sin*sin = cos
}	

// Yderive

/*! @brief Compute \f$ B \leftarrow D_y(A) \f$.
 *
 *	@note The basis functions are Fourier along y direction. So parity does not matter. 
 *
 *  @return \f$ B \leftarrow D_y(A) \f$.
 */
void IncVF::Yderiv_RSprod_VT(Array<complx,3> A, Array<complx,3> B)
{
	Yderiv(basis_type, N, A, B, kfactor, 0);
}


// Zderive	

/*! @brief Compute \f$ B \leftarrow D_z(A) \f$.
 *
 *	@note The basis functions are Fourier along z direction. So parity does not matter. 
 *
 *  @return \f$ B \leftarrow D_z(A) \f$.
 */
void IncVF::Zderiv_RSprod_VT(Array<complx,3> A, Array<complx,3> B)
{
	Zderiv(basis_type, N, A, B, kfactor, 0);
}

//****************************  End of Compute_derivative.cc **********************************






