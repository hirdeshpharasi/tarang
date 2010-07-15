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

/*! \file  compute_force_decay.cc
 * 
 * @brief  F = 0.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "../IncFluid.h"

//*********************************************************************************************


void IncFluid::Compute_force_decay()
{
	
	*Force1 = 0.0; 
	*Force2 = 0.0; 
	*Force3 = 0.0;
	
	if (alias_switch == "DEALIAS")		Dealias_force();	
}

//*********************************************************************************************

void IncFluid::Compute_force_decay(IncSF& T)
{
	
	Compute_force_decay();
		
	*T.Force = 0.0;
	
	if (alias_switch == "DEALIAS")		Dealias_force(T);	
}

//*********************************************************************************************


void IncFluid::Compute_force_decay(IncVF& W)
{
	
	Compute_force_decay();

	*W.Force1 = 0.0; 
	*W.Force2 = 0.0; 
	*W.Force3 = 0.0;
	
	if (alias_switch == "DEALIAS")		Dealias_force(W);
	
}


//*********************************************************************************************

void IncFluid::Compute_force_decay(IncVF& W, IncSF &T)
{
	Compute_force_decay(W);
	
	*T.Force = 0.0;
	
	if (alias_switch == "DEALIAS")		Dealias_force(W, T);
}

//************************ End of compute_force_decay.cc **************************************


