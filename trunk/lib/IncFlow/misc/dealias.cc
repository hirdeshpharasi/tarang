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

/*! \file  dealias.cc
 * 
 * @brief  Dealiase field, force, nlin
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */


#include "../IncVF.h"
#include "../IncSF.h"



/**********************************************************************************************

					Dealias IncVF

***********************************************************************************************/

void IncVF::Dealias()
{	
	CV_Dealias();
}

void IncVF::Dealias(IncVF& W)
{	
	CV_Dealias();
	W.CV_Dealias();
}

void IncVF::Dealias(IncSF& T)
{	
	CV_Dealias();
	T.CS_Dealias();
}

void IncVF::Dealias(IncVF& W, IncSF& T)
{	
	CV_Dealias();
	W.CV_Dealias();
}
/**********************************************************************************************

					Dealias forces

***********************************************************************************************/

void IncVF::Dealias_force()
{	

	Dealias_array(basis_type, N, *Force1);
	Dealias_array(basis_type, N, *Force2);
	Dealias_array(basis_type, N, *Force3);
}

void IncVF::Dealias_force(IncVF& W)
{	
	Dealias_force();
	W.Dealias_force();
}

void IncVF::Dealias_force(IncSF& T)
{	
	Dealias_force();
	
	Dealias_array(basis_type, N, *T.Force);
}


void IncVF::Dealias_force(IncVF& W, IncSF& T)
{	
	Dealias_force(T);
	W.Dealias_force();
}


/**********************************************************************************************

					Dealias nlin

***********************************************************************************************/

void IncVF::Dealias_nlin()
{	

	Dealias_array(basis_type, N, *nlin1);
	Dealias_array(basis_type, N, *nlin2);
	Dealias_array(basis_type, N, *nlin3);
}

void IncVF::Dealias_nlin(IncVF& W)
{	
	Dealias_nlin();
	W.Dealias_nlin();
}

void IncVF::Dealias_nlin(IncSF& T)
{	
	Dealias_nlin();
		
	Dealias_array(basis_type, N, *T.nlin);
}


void IncVF::Dealias_nlin(IncVF& W, IncSF& T)
{	
	Dealias_nlin(T);
	W.Dealias_nlin();
}


//*********************************  End of dealias.cc  ***************************************



	
