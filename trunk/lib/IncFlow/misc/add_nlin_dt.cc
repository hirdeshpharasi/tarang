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

/*! \file  add_nlin_dt.cc
 * 
 * @brief  V = V + nlin*dt;  F = F + nlin*dt
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

					V = V + nlin*dt;  F = F + nlin*dt

***********************************************************************************************/

void IncVF::Add_nlin_dt(DP dt)
{	
	*V1 = *V1 + complex<DP>(dt,0)*(*nlin1);           // V time-advance
	*V2 = *V2 + complex<DP>(dt,0)*(*nlin2);
	*V3 = *V3 + complex<DP>(dt,0)*(*nlin3);
}
 
//
//

void IncSF::Add_nlin_dt(DP dt)
{

	*F = *F + complex<DP>(dt,0)*(*nlin);
}

/**********************************************************************************************

					for CVF: W = W + factor*nlin; *T.F = *T.F + nlin*dt

***********************************************************************************************/

void IncVF::Add_nlin_to_field(CVF& W, DP factor)
{	
	*W.V1 = *W.V1 + complex<DP>(factor,0) * (*nlin1);         
	*W.V2 = *W.V2 + complex<DP>(factor,0) * (*nlin2);
	*W.V3 = *W.V3 + complex<DP>(factor,0) * (*nlin3); 
}
 
//
//

void IncSF::Add_nlin_to_field(CSF& T, DP factor)
{  
							  
	*T.F = *T.F + complex<DP>(factor,0) * (*nlin);
}


/**********************************************************************************************

					for PlainCVF: W = W + factor*nlin; *T.F = *T.F + nlin*dt

***********************************************************************************************/

void IncVF::Add_nlin_to_field(PlainCVF& W, DP factor)
{	
	*W.V1 = *W.V1 + complex<DP>(factor,0) * (*nlin1);         
	*W.V2 = *W.V2 + complex<DP>(factor,0) * (*nlin2);
	*W.V3 = *W.V3 + complex<DP>(factor,0) * (*nlin3); 
}
 
//
//

void IncSF::Add_nlin_to_field(PlainCSF& T, DP factor)
{  
							  
	*T.F = *T.F + complex<DP>(factor,0) * (*nlin);
}


//*****************************  End of add_nlin_dt.cc  ***************************************



	

