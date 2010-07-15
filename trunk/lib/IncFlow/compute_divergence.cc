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
/*! \file compute_divergence.h 
 * 
 * @brief Functions for computing divergence of nlin and V field.
 *
 * @note  Divergence is stored in CVF of class NLIN.  We use the same class to store
 *			pressure also.  Hence divergence function should be called in the beginning 
 *			or end of loop when *F is free.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date 30 August 2008
 * 
 * @bug  No known bugs
 */ 

//*********************************************************************************************


#include "IncVF.h"


/** @brief Compute divergence of nlin \f$ \vec{N} \f$.
 *
 *	@note	Divergence is stored in CVF of class NLIN.  We use the same class to store
 *			pressure also.  Hence divergence function should be called in the beginning 
 *			or end of loop when *F is free.
 *
 *  @return  \f$ *F = \mathcal{F}(D_i N_i) \f$. 
 */
 void IncVF::Compute_divergence_nlin()
{

	Xderiv(basis_type, N, *nlin1, *F, kfactor, 1);		
	
	Yderiv(basis_type, N, *nlin2, *VF_temp, kfactor, 0);
	*F = *F + *VF_temp;
	
	Zderiv(basis_type, N, *nlin3, *VF_temp, kfactor, 0);
	*F = *F + *VF_temp;
}

//********************************************************************************************* 




/** @brief Compute divergence of VF \f$ \vec{V} \f$.
 *
 *	@note	Divergence is stored in CVF of class NLIN.  We use the same class to store
 *			pressure also.  Hence divergence function should be called in the beginning 
 *			or end of loop when *F is free.
 *
 *  @return  \f$ *F = \mathcal{F}(D_i V_i) \f$. 
 */
void IncVF::Compute_divergence_field()
{

	Xderiv(basis_type, N, *V1, *F, kfactor, 1);			// factor =1 (sin)
	
	Yderiv(basis_type, N, *V2, *VF_temp, kfactor, 0);
	*F = *F + *VF_temp;
	
	Zderiv(basis_type, N, *V3, *VF_temp, kfactor, 0);
	*F = *F + *VF_temp;
	
	DP MYEPS2 = 10E-10;
	
	for (int i1=0; i1<local_N1; i1++)
		for (int i2=0; i2<N[2]; i2++)
			for (int i3=0; i3<=N[3]/2; i3++)
			{
				if (abs((*F)(i1,i2,i3)) > MYEPS2)
					cout << "div: proc_id, l1, l2, l3, div " << my_id << " " << 
						i1 << "  " << i2 << " " << i3 << " " << abs((*F)(i1,i2,i3)) << endl;
			}
}

//*******************************  End of compute_divergence.cc *******************************



