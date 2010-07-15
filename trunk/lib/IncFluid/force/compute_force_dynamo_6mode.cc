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

/*! \file  compute_force_dynamo_6mode.cc
 * 
 * @brief  Set up force for 6 mode dynamo model of Verma et al.
 *
 *  @note is_force_field_modes_read is used to read only once.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug No known bugs
 */

#include "../IncFluid.h"

//*********************************************************************************************


void IncFluid::Compute_force_DYNAMO_SIX_MODE(IncVF& W)
{

	if (is_force_field_modes_read == 0)
	{
		int k0 = ((int) (*force_field_para)(1));
		DP amp101 = (*force_field_para)(2);
		DP amp011 = (*force_field_para)(3);
		DP amp112 = (*force_field_para)(4);
		DP h = (*force_field_para)(5);
		
		Setup_SIX_MODE_force_field(k0, amp101, amp011, amp112, h); 
		// force only u field
		
		if (my_id == master_id)
		{
			cout	<< "FORCE 101  " 
					<< (*Force1)(1,0,1) << "  " 
					<< (*Force2)(1,0,1) << "  "
					<< (*Force3)(1,0,1) << endl;
					
			cout	<< "FORCE 011  " 
					<< (*Force1)(0,1,1) << "  " 
					<< (*Force2)(0,1,1) << "  " 
					<< (*Force3)(0,1,1) << endl;
					
			cout	<< "FORCE 112  " 
					<< (*Force1)(1,1,2) << "  " 
					<< (*Force2)(1,1,2) << "  " 
					<< (*Force3)(1,1,2) << endl;
		}			
		
		is_force_field_modes_read = 1;
	}
			
}



//****************************** End of compute_force_dynamo_6mode.cc *************************


