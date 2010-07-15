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

/*! \file  init_cond_dynamo_6mode.cc
 * 
 * @brief Initial conditions for dynamo 6 model of Verma et al.
 *
 * @note The parameters are read from parameter file. 
 *		  Applicable for FOUR basis; not for SCFT basis.
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug No known bugs
 */

#include "../IncFluid.h"

//*********************************************************************************************

void IncFluid::Init_cond_DYNAMO_SIX_MODE(IncVF& W)
{

	if (basis_type == "FOUR")
	{
		int k0 = ((int) (*init_cond_para)(1));
		DP amp101 = (*init_cond_para)(2);
		DP amp011 = (*init_cond_para)(3);
		DP amp112 = (*init_cond_para)(4);
		DP h = (*init_cond_para)(5);
		
		DP ampW101 = (*init_cond_para)(6);
		DP ampW011 = (*init_cond_para)(7);
		DP ampW112 = (*init_cond_para)(8);
		DP hW = (*init_cond_para)(9);
		
		Setup_SIX_MODE_field(k0, amp101, amp011, amp112, h);
		W.Setup_SIX_MODE_field(k0, ampW101, ampW011, ampW112, hW);
	}
	
	else if (basis_type == "SCFT")
	{
		cout << "Setup_SIX_MODE_field Not allowed in SCFT basis_type " << endl;
		exit(1);
	}
}

//******************************* End of Init_cond_dynamo_6mode.cc ****************************


