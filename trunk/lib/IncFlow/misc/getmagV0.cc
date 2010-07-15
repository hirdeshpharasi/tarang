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


/*! \file  Get_mag_V0.cc
 * 
 * @brief  Returns the magnitude of the mean field.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */


#include "../IncVF.h"


//*********************************************************************************************	

DP	IncVF::Get_mag_V0()
{
	DP modV0;
	
	if (my_id == master_id)
		modV0 = sqrt( pow2(real((*V1)(0,0,0)))  +  pow2(real((*V2)(0,0,0)))  
									+  pow2(real((*V3)(0,0,0))) ) ;
									
	MPI_Bcast( &modV0, 1, MPI_DOUBLE, master_id, MPI_COMM_WORLD);	
	
	return modV0;														
}



//*****************************  End of copy_field.cc  ****************************************


	
