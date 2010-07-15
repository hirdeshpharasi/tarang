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
/*! \file  compute_force_RB.cc
 * 
 * @brief Force for RB convection
 *
 *		For Pr_large, u_small:		 F1 = Ra*Pr*T(k);  
 *		For Pr_large, u_large:		 F1 = Ra*T(k);    	
 *		For Pr_small or 0, u_small:  F1 = Ra*T(k);	   
 *		For Pr_small or 0, u_large:  F1 = Pr*T(k)
 *
 *		For Pr_large				Ftheta = V1;    	
 *		For Pr_small				Ftheta = V1/Pr;	   
 *		For Pr = 0					Ftheta = 0.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */
 
#include "../IncFluid.h"


//*********************************************************************************************

void IncFluid::Compute_force_RB(IncSF& T, DP Ra, DP Pr, string Pr_switch, string RB_Uscaling)
{

// For the velocity field
	if (Pr_switch == "PRLARGE") 
	{
	
		if (RB_Uscaling == "USMALL") 		
			*Force1 =  (Ra*Pr)*(*T.F);			// (u.grad)u-Ra*Pr*theta	
			
		else if (RB_Uscaling  == "ULARGE") 		
			 *Force1 =  (*T.F);					// (u.grad)u-theta			
	}
	
	else if ((Pr_switch == "PRSMALL") || (Pr_switch == "PRZERO")) 
	{
	
		if (RB_Uscaling == "USMALL") 
			*Force1 = (Ra)*(*T.F);				// (u.grad)u-theta
			
		else if (RB_Uscaling == "ULARGE") 		
			*Force1 = (Pr)*(*T.F);				// (u.grad)u-theta
	}	


	*Force2 = 0.0; 
	*Force3 = 0.0;

	
// For the temperature field

	if (Pr_switch == "PRZERO") 
		*T.Force = 0.0;  
		// Do nothing; Nonlinear term (u.grad)T does not exist- see single-time step
	
	else if (Pr_switch == "PRLARGE")	
		*T.Force =   *V1;										// F(T) = ux(k)	
	
	else if (Pr_switch == "PRSMALL") 	
		*T.Force =   complex<DP>(1/Pr, 0)*((*V1)); ;			//F(T) =  ux(k)/Pr


	if (alias_switch == "DEALIAS")		Dealias_force(T);	
	
			
}



//*********************************************************************************************


void IncFluid::Compute_force_RB
(
	IncVF& W, IncSF& T, 
	DP Ra, DP Pr, 
	string Pr_switch, string RB_Uscaling
)
{
	
	Compute_force_RB(T, Ra, Pr, Pr_switch, RB_Uscaling);

	*W.Force1 = 0.0; 
	*W.Force2 = 0.0; 
	*W.Force3 = 0.0;	

	if (alias_switch == "DEALIAS")		W.Dealias_force();
}


//************************ End of compute_force_RB.cc *****************************************

