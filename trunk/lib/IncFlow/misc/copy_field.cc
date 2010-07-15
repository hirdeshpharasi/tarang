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


/*! \file  copy_field.cc
 * 
 * @brief  Copy field
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
	For PlainCVF & PlainCSF
					Copy_field_to(W):	  W.V <- V
					Copy_field_from(W):   V <-W. V

***********************************************************************************************/

void IncVF::Copy_field_to(PlainCVF& W)
{	

	*W.V1 = *V1;      
	*W.V2 = *V2;		
	*W.V3 = *V3;
}

//
//
void IncSF::Copy_field_to(PlainCSF& T)
{
							  
	*T.F = *F;
}  

//
//

void IncVF::Copy_field_from(PlainCVF& W)
{	

	*V1 = *W.V1;      
	*V2 = *W.V2;		
	*V3 = *W.V3;
}

//
//
void IncSF::Copy_field_from(PlainCSF& T)
{
							  
	*F = *T.F;
}  



/**********************************************************************************************
		For CVF & CSF
					Copy_field_to(W):	  W.V <- V
					Copy_field_from(W):   V <-W. V

***********************************************************************************************/

void IncVF::Copy_field_to(CVF& W)
{	

	*W.V1 = *V1;      
	*W.V2 = *V2;		
	*W.V3 = *V3;
}

//
//
void IncSF::Copy_field_to(CSF& T)
{
							  
	*T.F = *F;
}  

//
//

void IncVF::Copy_field_from(CVF& W)
{	

	*V1 = *W.V1;      
	*V2 = *W.V2;		
	*V3 = *W.V3;
}

//
//
void IncSF::Copy_field_from(CSF& T)
{
							  
	*F = *T.F;
}  



//*****************************  End of copy_field.cc  ****************************************

