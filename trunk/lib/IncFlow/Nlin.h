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

/*! \file  Nlin.h
 * 
 * @brief contains nonlinear arrays \f$ \vec{N} = \vec{u}\cdot \nabla \vec{u} \f$.
 *	Inherits a scalar field for pressure
 * 
 * @sa IncVF.h
 * @sa IncSF.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug   No known bugs
 */

//*********************************************************************************************

#ifndef _NLIN_H
#define _NLIN_H

#include "../fields/fields.h"


// nonlinear field; Inherits a scalar
class NLIN : public CSF
{
 public:

	//! For storing T[U.grad U]
	Array<complx,3> *nlin1;   
	Array<complx,3> *nlin2;   
	Array<complx,3> *nlin3;  
											
	int		Nn[4];
	string	NLIN_basis_type;
	DP		NLIN_kfactor[4];

 
//*********************************************************************************************
 
public:

	NLIN
	(
		int *NN, 
		string string_switches[],
		Array<int,1> switches, 
		DP *kfactor, 
		Array<int,1> misc_output_para
	);  
  
  
	/*! @brief	After the execution \f$ \vec{N} = \mathcal{F}(D_i V_i^2) \f$, where
	 * \f$ D_i \f$ is the derivative along the i-th direction.
	 *
	 * This is the diagonal term of the \f$ \nabla \vec{V} \vec{V} \f$.
	 *  Note that \f$ V_i^2 \f$ have even parity in SCFT basis.
	 *
	 * @note Before entering this function \f$ \vec{N} = V_i^2 \f$,
	 */
	void NLIN_diag_Forward_transform_derivative(Array<complx,3> temp_r);									
};

#endif

//******************************  End Class decl of NLIN   ************************************


