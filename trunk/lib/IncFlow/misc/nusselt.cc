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
/*! \file  nusselt.cc
 * 
 * @brief Computes nusselt number
 *
 * @author  M. K. Verma
 * @version 4.0  MPI
 * @date Feb 2009
 *
 * @bug  No known bugs
 */

#include "../IncVF.h"
#include "../IncSF.h"


//*********************************************************************************************

/** @brief Computes Nusselt number
 * 
 *  @param IncSF T
 *  @param Ra  Rayleigh number
 *  @param Pr  Prandtl number
 *  @param  string Pr_switch, string RB_Uscaling
 *
 * @return \f$ Nu = 1 + \sum u_z T \f$ 
 * @return For Pr=0, Nu = 1, but the function returns \f$ Nu = \sum u_z T \f$
 */
DP IncVF::Get_Nusselt_no(IncSF& T)
{

	// Actually Nu = 1 for Pr=0.  Here we just report the product.
	if (globalvar_Pr_switch == "PRZERO") 
		return ( 2* Get_total_energy(basis_type, alias_switch, N, *V1, *T.F) );
	
	
	else if (globalvar_Pr_switch == "PRLARGE") 
	{
		if (globalvar_RB_Uscaling == "USMALL") 
			return ( 1 + 2* Get_total_energy(basis_type, alias_switch, N, *V1, *T.F) );
			
		else if (globalvar_RB_Uscaling == "ULARGE") 
			return ( 1 + 2*sqrt(globalvar_Ra*globalvar_Pr)
							* Get_total_energy(basis_type, alias_switch,N, *V1, *T.F) );
	}
	
	else if (globalvar_Pr_switch == "PRSMALL") 
	{
		if (globalvar_RB_Uscaling == "USMALL") 
			return ( 1 + 2*globalvar_Pr*globalvar_Pr
					* Get_total_energy(basis_type, alias_switch, N, *V1, *T.F) );
			
		else if (globalvar_RB_Uscaling == "ULARGE") 
			return ( 1 + 2*globalvar_Pr*sqrt(globalvar_Ra*globalvar_Pr)
					* Get_total_energy(basis_type, alias_switch, N, *V1, *T.F) );
	}
	
	else if (globalvar_Pr_switch == "PRINFTY") 
		return ( 1 + 2* Get_total_energy(basis_type, alias_switch, N, *V1, *T.F) );

	
	else
		return 0;	// forw -Wall

}


//***********************************  End of nusselt.cc **************************************



