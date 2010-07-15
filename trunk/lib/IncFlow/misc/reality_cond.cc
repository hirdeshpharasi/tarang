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

/*! \file  realitiy_cond.cc
 * 
 * @brief Functions to satisfy the reality conditions for the fields and the force. 
 *			Fill -K vectors in ky=0 & N[2]/2 line in 2D, and kz=0 & N[3]/2 planes in 3D.
 *
 * @sa void IncVF::Satisfy_reality_condition_field()
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug	  No known bugs
 */

#include "../IncVF.h"
#include "../IncSF.h"

//*********************************************************************************************	

void IncVF::Satisfy_reality_condition_field()
{
	Satisfy_reality_array(basis_type, N, *V1);
	Satisfy_reality_array(basis_type, N, *V2);
	Satisfy_reality_array(basis_type, N, *V3);
}

//
//

void IncVF::Satisfy_reality_condition_field(IncSF& T)
{

	Satisfy_reality_condition_field();				// For V field
	
	Satisfy_reality_array(basis_type, N, *T.F);
}

//
//

void IncVF::Satisfy_reality_condition_field(IncVF& W)
{
	Satisfy_reality_condition_field();
	W.Satisfy_reality_condition_field();
}

//
//


void IncVF::Satisfy_reality_condition_field(IncVF& W, IncSF& T)
{

	Satisfy_reality_condition_field();				// For V field
	W.Satisfy_reality_condition_field();
	
	Satisfy_reality_array(basis_type, N, *T.F);	
}



//*********************************************************************************************	

void IncVF::Satisfy_reality_condition_force_field()
{
	Satisfy_reality_array(basis_type, N, *Force1);
	Satisfy_reality_array(basis_type, N, *Force2);
	Satisfy_reality_array(basis_type, N, *Force3);
}

//
//

void IncVF::Satisfy_reality_condition_force_field(IncSF& T)
{
	
	Satisfy_reality_condition_force_field();				// For V field
	
	Satisfy_reality_array(basis_type, N, *T.Force);
}

//
//

void IncVF::Satisfy_reality_condition_force_field(IncVF& W)
{
	Satisfy_reality_condition_force_field();
	W.Satisfy_reality_condition_force_field();
}

//
//


void IncVF::Satisfy_reality_condition_force_field(IncVF& W, IncSF& T)
{
	
	Satisfy_reality_condition_force_field();				// For V field
	W.Satisfy_reality_condition_force_field();
	
	Satisfy_reality_array(basis_type, N, *T.Force);	
}

//*********************************************************************************************	

void IncVF::Satisfy_reality_condition_nlin()
{
	Satisfy_reality_array(basis_type, N, *nlin1);
	Satisfy_reality_array(basis_type, N, *nlin2);
	Satisfy_reality_array(basis_type, N, *nlin3);
}

//
//

void IncVF::Satisfy_reality_condition_nlin(IncSF& T)
{
	
	Satisfy_reality_condition_nlin();				// For V field
	
	Satisfy_reality_array(basis_type, N, *T.nlin);
}

//
//

void IncVF::Satisfy_reality_condition_nlin(IncVF& W)
{
	Satisfy_reality_condition_nlin();
	W.Satisfy_reality_condition_nlin();
}

//
//


void IncVF::Satisfy_reality_condition_nlin(IncVF& W, IncSF& T)
{
	
	Satisfy_reality_condition_nlin();				// For V field
	W.Satisfy_reality_condition_nlin();
	
	Satisfy_reality_array(basis_type, N, *T.nlin);	
}

/**********************************************************************************************	

		Satisfy reality condition for the fields + force field + nlin
		In the egdes..
		
***********************************************************************************************/

void IncVF::Satisfy_reality_condition()
{

	Satisfy_reality_array(basis_type, N, *V1);
	Satisfy_reality_array(basis_type, N, *V2);
	Satisfy_reality_array(basis_type, N, *V3);
	
	Satisfy_reality_array(basis_type, N, *Force1);
	Satisfy_reality_array(basis_type, N, *Force2);
	Satisfy_reality_array(basis_type, N, *Force3);
	
	Satisfy_reality_array(basis_type, N, *nlin1);
	Satisfy_reality_array(basis_type, N, *nlin2);
	Satisfy_reality_array(basis_type, N, *nlin3);
}

//
//

void IncVF::Satisfy_reality_condition(IncSF& T)
{

	Satisfy_reality_condition();				// For V field
	
	// for scalar fn
	Satisfy_reality_array(basis_type, N, *T.F);
	
	Satisfy_reality_array(basis_type, N, *T.Force);
		
	Satisfy_reality_array(basis_type, N, *T.nlin);	

}

//
//

void IncVF::Satisfy_reality_condition(IncVF& W)
{
	Satisfy_reality_condition();
	W.Satisfy_reality_condition();
}

//
//


void IncVF::Satisfy_reality_condition(IncVF& W, IncSF& T)
{

	Satisfy_reality_condition();				// For V field
	W.Satisfy_reality_condition();
	
	Satisfy_reality_array(basis_type, N, *T.F);
	
	Satisfy_reality_array(basis_type, N, *T.Force);
		
	Satisfy_reality_array(basis_type, N, *T.nlin);
	
}



//**************************** End of reality_cond.cc *****************************************	



