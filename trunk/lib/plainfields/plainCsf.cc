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


/*! \file  PlainCsf.cc
 * 
 * @brief  Class declaration of Csf, a Complex Scalar Field 
 * @sa PlainCsf.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug	 No known bugs
 */

#include "plainCsf.h"


//*********************************************************************************************

PlainCSF::PlainCSF(int *NN)
{

	for (int i=1; i<=3; i++) 
		Ncs_plain[i] = NN[i];  
	
 
	F=new Array<complx,3>(local_N1, Ncs_plain[2], Ncs_plain[3]/2+1);  
	*F=0.0;   // initialize
}

//*********************************************************************************************

void PlainCSF::PlainCS_Initialize()
{
	
	*F = 0.0;
}

//******************************** End of CSF.cc **********************************************


