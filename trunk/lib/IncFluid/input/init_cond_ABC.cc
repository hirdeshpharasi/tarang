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


/*! \file  init_cond_ABC.cc
 * 
 * @brief Initial conditions as ABC flow.
 *
 * @note The parameters are read from parameter file. 
 *		  Applicable for FOUR basis; not for SCFT basis.
 *
 * @note \f$ V_x = amp (B \cos(k_0 y) + C \sin(k_0 z)) \f$
 *		 \f$ V_y = amp (A \sin(k_0 x) + C \cos(k_0 z)) \f$
 *		 \f$ V_z = amp (A \cos(k_0 x) + C \cos(k_0 y)) \f$
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */
 

#include "../IncFluid.h"


//*********************************************************************************************

void IncFluid::Init_cond_ABC()
{

	if (basis_type == "FOUR")
	{
		int k0 = ((int) (*init_cond_int_para)(1));
		DP A = (*init_cond_int_para)(2);
		DP B = (*init_cond_int_para)(3);
		DP C = (*init_cond_int_para)(4);
		DP amp = (*init_cond_int_para)(5);
		
		Setup_ABC_field(k0, amp, A, B, C);
	}	
	
	else if (basis_type == "SCFT")
	{
		cout << "ABC field not allowed for SCFT " << endl;
		exit(1);
	}
	
}

//*********************************************************************************************
// Scalar 
//

void IncFluid::Init_cond_ABC(IncSF& T)
{
	
	if (basis_type == "FOUR")
	{
		Init_cond_ABC();
		
		*T.F = 0.0;
	}
	
	else if (basis_type == "SCFT")
	{
		cout << "ABC field not allowed for SCFT " << endl;
		exit(1);
	}
	
}


//*********************************************************************************************
// Vector
//

void IncFluid::Init_cond_ABC(IncVF& W)
{

	if (basis_type == "FOUR")
	{
		int k0 = ((int) (*init_cond_int_para)(1));
		DP A = (*init_cond_int_para)(2);
		DP B = (*init_cond_int_para)(3);
		DP C = (*init_cond_int_para)(4);
		DP amp = (*init_cond_int_para)(5);
		DP ampW = (*init_cond_int_para)(6);
		
		Setup_ABC_field(k0, amp, A, B, C);
		W.Setup_ABC_field(k0, ampW, A, B, C);
	}
	
	else if (basis_type == "SCFT")
	{
		cout << "ABC field not allowed for SCFT " << endl;
		exit(1);
	}
		

}

//*********************************************************************************************
//	Vector+scalar
//

void IncFluid::Init_cond_ABC(IncVF& W, IncSF& T)
{

	if (basis_type == "FOUR")
	{
		Init_cond_Taylor_Green(W);
		
		*T.F = 0.0;
	}
	
	else if (basis_type == "SCFT")
	{
		cout << "ABC field not allowed for SCFT " << endl;
		exit(1);
	}

}



//******************************** End of Init_cond_ABC.cc ************************************

