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

/*! \file elsasser.cc
 * 
 * CONVERTS FIELD VARS AND NLINS FROM U_B TO ELSASSER VARS AND VICE VERSA
 *		
 *	void UB_to_Elsasser_field(IncVF& W) -- U(Zp) = U + W; W(Zm) = U - W;
 *	void Elsasser_to_UB_field(IncVF& W) -- U = (U(Zp) + W(Zm))/2; W =  (U(Zp) - W(Zm))/2; 
 *	
 *	Same as above for force.
 *	
 *	void UB_to_Elsasser_nlin(IncVF& W) 
 *		-- U.nlin(Zp) = U.nlin + W.nlin; W(Zm) = U.nlin - W.nlin;
 *	void Elsasser_to_UB(IncVF& W) 
 *		-- U.nlin = (U(Zp).nlin + W(Zm).nlin)/2; W =  (U.nlin(Zp) - W.nlin(Zm))/2; 
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug   No known bugs
 */

#include "../IncVF.h"
#include "../IncSF.h"

//*********************************************************************************************

void IncVF::UB_to_Elsasser_field(IncVF& W)
{
	*VF_temp = *V1;   *V1 = *V1 + *W.V1;   *W.V1 = *VF_temp - *W.V1; 
	
	*VF_temp = *V2;   *V2 = *V2 + *W.V2;   *W.V2 = *VF_temp - *W.V2;
	
	*VF_temp = *V3;   *V3 = *V3 + *W.V3;   *W.V3 = *VF_temp - *W.V3; 
}

void IncVF::Elsasser_to_UB_field(IncVF& W)
{
	*VF_temp = *V1;   
	*V1 =  complex<double>(0.5,0)*(*V1 + *W.V1);   
	*W.V1 =  complex<double>(0.5,0)*(*VF_temp - *W.V1);
	 
	*VF_temp = *V2;   
	*V2 =  complex<double>(0.5,0)*(*V2 + *W.V2);   
	*W.V2 =  complex<double>(0.5,0)*(*VF_temp - *W.V2); 
	
	*VF_temp = *V3;   
	*V3 =  complex<double>(0.5,0)*(*V3 + *W.V3);   
	*W.V3 =  complex<double>(0.5,0)*(*VF_temp - *W.V3);
}

// Force field
//

void IncVF::UB_to_Elsasser_force(IncVF& W)
{
	*VF_temp = *Force1;   *Force1 = *Force1 + *W.Force1;   *W.Force1 = *VF_temp - *W.Force1; 
	
	*VF_temp = *Force2;   *Force2 = *Force2 + *W.Force2;   *W.Force2 = *VF_temp - *W.Force2;
	
	*VF_temp = *Force3;   *Force3 = *Force3 + *W.Force3;   *W.Force3 = *VF_temp - *W.Force3; 
}

void IncVF::Elsasser_to_UB_force(IncVF& W)
{
	*VF_temp = *Force1;   
	*Force1 =  complex<double>(0.5,0)*(*Force1 + *W.Force1);   
	*W.Force1 =  complex<double>(0.5,0)*(*VF_temp - *W.Force1); 
	
	*VF_temp = *Force2;   
	*Force2 =  complex<double>(0.5,0)*(*Force2 + *W.Force2);   
	*W.Force2 =  complex<double>(0.5,0)*(*VF_temp - *W.Force2); 
	
	*VF_temp = *Force3;   
	*Force3 =  complex<double>(0.5,0)*(*Force3 + *W.Force3);   
	*W.Force3 =  complex<double>(0.5,0)*(*VF_temp - *W.Force3);
}

//
// nlin
//

void IncVF::UB_to_Elsasser_nlin(IncVF& W)
{
	*VF_temp = *nlin1;   *nlin1 = *nlin1 + *W.nlin1;   *W.nlin1 = *VF_temp - *W.nlin1; 
	
	*VF_temp = *nlin2;   *nlin2 = *nlin2 + *W.nlin2;   *W.nlin2 = *VF_temp - *W.nlin2; 
	
	*VF_temp = *nlin3;   *nlin3 = *nlin3 + *W.nlin3;   *W.nlin3 = *VF_temp - *W.nlin3;
}

void IncVF::Elsasser_to_UB_nlin(IncVF& W)
{
	*VF_temp = *nlin1;   
	*nlin1 = complex<double>(0.5,0)*(*nlin1 + *W.nlin1);   
	*W.nlin1 = complex<double>(0.5,0)*(*VF_temp - *W.nlin1); 
	
	*VF_temp = *nlin2;   
	*nlin2 = complex<double>(0.5,0)*(*nlin2 + *W.nlin2);   
	*W.nlin2 = complex<double>(0.5,0)*(*VF_temp - *W.nlin2); 
	
	*VF_temp = *nlin3;   
	*nlin3 = complex<double>(0.5,0)*(*nlin3 + *W.nlin3);   
	*W.nlin3 = complex<double>(0.5,0)*(*VF_temp - *W.nlin3); 
}

//******************************  End of Elasser.cc  ******************************************

	
