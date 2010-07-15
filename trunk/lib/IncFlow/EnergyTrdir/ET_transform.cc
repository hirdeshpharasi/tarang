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

/*! \file  ET_transform.cc
 * 
 * @brief  Forward and inverse transforms for energy transfer functions
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "../IncVF.h"
					

//*********************************************************************************************

/*! @brief Inplace inverse transform of Giver field G.
 *
 * @note temp_r is temporary array.
 *
 * @return Inplace inverse transform of Giver field G.
 */
void IncVF::ET_Inverse_transform_G(Array<complx,3> temp_r)
{					
	Inverse_transform_array(basis_type, N, *G1, temp_r, 1);			// ISFT
	Inverse_transform_array(basis_type, N, *G2, temp_r, 0);			// ICFT
	Inverse_transform_array(basis_type, N, *G3, temp_r, 0);			// ICFT
}		


//*********************************************************************************************

/*! @brief inverse transform of Giver field G under transpose order
 *
 * @return Inplace inverse transform of Giver field G.
 */
 
void IncVF::ET_Inverse_transform_G_transpose_order(Array<complx,3> temp)
{	
#ifdef TRANSPOSE	
	temp = *G1;
	Inverse_transform_array_transpose_order(basis_type, N, temp, *G1, 1);	// ISFT
	
	temp = *G2;
	Inverse_transform_array_transpose_order(basis_type, N, temp, *G2, 0);	// ICFT
	
	temp = *G3;
	Inverse_transform_array_transpose_order(basis_type, N, temp, *G3, 0);	// ICFT
#endif	
}		


//*********************************************************************************************

/*! @brief Inplace forward transform of Giver field G.
 *
 * @note temp_r is temporary array.
 *
 * @return Inplace inverse transform of Giver field G.
 */
void IncVF::ET_Forward_transform_RSprod_G(Array<complx,3> temp_r)
{
	Forward_transform_array(basis_type, N, *G1, temp_r, 1);			// SFT
	Forward_transform_array(basis_type, N, *G2, temp_r, 0);			// CFT
	Forward_transform_array(basis_type, N, *G3, temp_r, 1);			// SFT
}


//*********************************************************************************************

/*! @brief forward transform of Giver field G under transpose order
 *
 * @return forward transform of Giver field G.
 */
void IncVF::ET_Forward_transform_RSprod_G_transpose_order(Array<complx,3> temp_r)
{
#ifdef TRANSPOSE
	temp_r = *G1;
	Forward_transform_array_transpose_order(basis_type, N, temp_r, *G1, 1);		// SFT
	
	temp_r = *G2;
	Forward_transform_array_transpose_order(basis_type, N, temp_r, *G2, 0);		// CFT
	
	temp_r = *G3;
	Forward_transform_array_transpose_order(basis_type, N, temp_r, *G3, 1);		// SFT
#endif	
}

//********************************  End of ET_transform.cc ************************************





