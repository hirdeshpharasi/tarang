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

/*! \file  Nlin.cc
 * 
 * @sa Nlin.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @Bug  No known bugs
 */
 
#include "Nlin.h"
// #include "../fields/fields.h"
#include <complex>
					
					
//*********************************************************************************************

NLIN::NLIN
(
	int *NN, 
	string string_switches[], 
	Array<int,1> switches, 
	DP *kfactor, 
	Array<int,1> misc_output_para
):		CSF(NN, string_switches, switches, kfactor, misc_output_para) 
{

	NLIN_basis_type = string_switches[1];
	
	for (int i=1; i<=3; i++) 
	{
		Nn[i]=NN[i];	
		NLIN_kfactor[i] = kfactor[i];
    }  
   
   
	nlin1 = new Array<complx,3>(local_N1, Nn[2], Nn[3]/2+1);  
	nlin2 = new Array<complx,3>(local_N1, Nn[2], Nn[3]/2+1);  
	nlin3 = new Array<complx,3>(local_N1, Nn[2], Nn[3]/2+1); 

	*nlin1 = 0.0; 
	*nlin3 = 0.0;  
	*nlin3 = 0.0;
}


//*********************************************************************************************

/*! @brief Compute \f$ N_i \leftarrow D_i \mathcal{F}(Vr_i * Wr_i) \f$ [no i sum].
 *
 *	For fluid, \f$ Wr_i = Vr_i \f$ 
 *
 *	Before entering the function <BR>
 *		\f$ N_i = \mathcal{F}(Vr_i * Wr_i) \f$ <BR>
 *
 *  @return \f$ N_i \leftarrow D_i \mathcal{F}(Vr_i * Wr_i) \f$.
 */
void NLIN::NLIN_diag_Forward_transform_derivative(Array<complx,3> temp_r)
{
	Forward_transform_array(NLIN_basis_type, Nn, *nlin1, temp_r, 0);
	Forward_transform_array(NLIN_basis_type, Nn, *nlin2, temp_r, 0);
	Forward_transform_array(NLIN_basis_type, Nn, *nlin3, temp_r, 0);
	
	Xderiv(NLIN_basis_type, Nn, *nlin1, *nlin1, NLIN_kfactor, 0);  
	Yderiv(NLIN_basis_type, Nn, *nlin2, *nlin2, NLIN_kfactor, 0);
	Zderiv(NLIN_basis_type, Nn, *nlin3, *nlin3, NLIN_kfactor, 0);	
}


//****************************** End Class decl of Nlin.c  ************************************




