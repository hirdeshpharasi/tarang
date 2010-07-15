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

/*! \file  add_force.cc
 * 
 * @brief Add force to nlin 
 *
 * @note	Before entering the function, nlin = V.grad(V). So nlin = nlin -Force.
 *			In Compute_rhs nlin = -nlin -grad(p).
 *
 * @sa void IncFluid::Add_force()
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */

#include "IncFluid.h"

//*********************************************************************************************

void IncFluid::Add_force()
{

	*nlin1 = *nlin1 - *Force1;	
	*nlin2 = *nlin2 - *Force2;
	*nlin3 = *nlin3 - *Force3;		
}

// Passive scalar or RB convection

void IncFluid::Add_force(IncSF& T)
{

	Add_force();

	*T.nlin = *T.nlin - *T.Force;						
}

// MHD

void IncFluid::Add_force(IncVF& W)
{

	Add_force();
	
	*W.nlin1 = *W.nlin1 - *W.Force1;	
	*W.nlin2 = *W.nlin2 - *W.Force2;
	*W.nlin3 = *W.nlin3 - *W.Force3;						
}

// MHD + Scalar

void IncFluid::Add_force(IncVF& W, IncSF& T)
{
	Add_force(W);
	
	*T.nlin = *T.nlin - *T.Force;						
}

//		RB Convection	//

void IncFluid::Add_force(IncSF& T, string Pr_switch)
{
	if (Pr_switch == "PRZERO")
		Add_force();
	else
		Add_force(T);
}		

void IncFluid::Add_force(IncVF& W, IncSF& T, string Pr_switch)
{
	if (Pr_switch == "PRZERO")
		Add_force(W);
	else
		Add_force(W, T);
}	


//******************************* End ofadd_force.cc ******************************************


