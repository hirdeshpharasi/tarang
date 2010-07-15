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

/*! \file  PlainCvf.cc
 * 
 * @brief  Class declaration of Cvf, a Complex Vector Field 
 *
 * @sa Cvf.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug  No known bugs
 */

#include "plainCvf.h"


/**********************************************************************************************

   Class constructor: Allocates for Vi; Initializes the arrays.
   
**********************************************************************************************/


PlainCVF::PlainCVF(int *NN)
{
	
	for (int i=1; i<=3; i++) 
		Ncv_plain[i]=NN[i];  


	V1 = new Array<complx,3>(local_N1, Ncv_plain[2], Ncv_plain[3]/2+1);  
	V2 = new Array<complx,3>(local_N1, Ncv_plain[2], Ncv_plain[3]/2+1);   
	V3 = new Array<complx,3>(local_N1, Ncv_plain[2], Ncv_plain[3]/2+1);   

	(*V1)=0.0; 
	(*V2)=0.0; 
	(*V3)=0.0;
}


/**********************************************************************************************

		 Initialize CVF

**********************************************************************************************/

void PlainCVF::PlainCV_Initialize()
{
	*V1 = 0.0;   
	*V2 = 0.0;   
	*V3 = 0.0;
}


//************************ END of PlainCVF class Definitions **********************************




