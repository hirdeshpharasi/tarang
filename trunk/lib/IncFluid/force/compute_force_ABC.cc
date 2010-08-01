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
/*! \file  compute_force_ABC.cc
 * 
 * @brief  Set up force using ABC model.
 *
 * @note The parameters are read from parameter file. 
 *		  Applicable for FOUR basis; not for SCFT basis.
 *
 * @note \f$ F_x = amp (B \cos(k_0 y) + C \sin(k_0 z)) \f$
 *		 \f$ F_y = amp (A \sin(k_0 x) + C \cos(k_0 z)) \f$
 *		 \f$ F_z = amp (A \cos(k_0 x) + C \cos(k_0 y)) \f$
 *
 *  @note is_force_field_para_read is used to read only once.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */


#include "../IncFluid.h"


//*********************************************************************************************

void IncFluid::Compute_force_ABC()
{
	
	if (is_force_field_para_read == 0)
	{
		if (basis_type == "FOUR")
		{
			int k0 = ((int) (*force_field_para)(1));
			DP A = (*force_field_para)(2);
			DP B = (*force_field_para)(3);
			DP C = (*force_field_para)(4);
			DP force_amp = (*force_field_para)(5);
			
			Setup_ABC_force_field(k0, force_amp, A, B, C);
		
			is_force_field_para_read = 1;
		}
		
		else if (basis_type == "SCFT")
		{
			cout << "ABC force field not allowed for SCFT " << endl;
			exit(1);
		}		
	}	
	
}


//*********************************************************************************************
//
// Scalar
//

void IncFluid::Compute_force_ABC(IncSF& T)
{

	if (is_force_field_para_read == 0)
	{
		if (basis_type == "FOUR")
		{
			*T.Force = 0.0;
			
			Compute_force_ABC();
			
			is_force_field_para_read = 1;
		}
		
		else if (basis_type == "SCFT")
		{
			cout << "ABC force field not allowed for SCFT " << endl;
			exit(1);
		}
				
	}

}


//*********************************************************************************************
//
// Force both V and W
//

void IncFluid::Compute_force_ABC(IncVF& W)
{	

	if (is_force_field_para_read == 0)
	{
		if (basis_type == "FOUR")
		{
			int k0 = ((int) (*force_field_para)(1));
			DP A = (*force_field_para)(2);
			DP B = (*force_field_para)(3);
			DP C = (*force_field_para)(4);
			DP force_amp = (*force_field_para)(5);
			DP forceW_amp = (*force_field_para)(6);
			
			Setup_ABC_force_field(k0, force_amp, A, B, C);
			W.Setup_ABC_force_field(k0, forceW_amp, A, B, C);
			
			is_force_field_para_read = 1;
		}
		
		else if (basis_type == "SCFT")
		{
			cout << "ABC force field not allowed for SCFT " << endl;
			exit(1);
		}		
	}	
	
}


//*********************************************************************************************
//
// Vector + scalar
//

void IncFluid::Compute_force_ABC(IncVF& W, IncSF& T)
{	

	if (is_force_field_para_read == 0)
	{
		if (basis_type == "FOUR")
		{
			*T.Force = 0.0;
			Compute_force_ABC(W);
			
			is_force_field_para_read = 1;
		}
		
		else if (basis_type == "SCFT")
		{
			cout << "ABC force field not allowed for SCFT " << endl;
			exit(1);
		}		
	}

}	

//*****************************   End of compute_force_ABC.cc *********************************



