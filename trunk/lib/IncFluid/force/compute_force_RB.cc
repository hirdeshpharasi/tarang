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

void IncFluid::Compute_force_RB(IncSF& T)
{

// For the velocity field
	if (globalvar_Pr_switch == "PRLARGE") {
		
		if (globalvar_RB_Uscaling == "USMALL") 		
			*Force1 =  (globalvar_Ra*globalvar_Pr)*(*T.F);			// (u.grad)u-Ra*Pr*theta	
		
		else if (globalvar_RB_Uscaling  == "ULARGE") 		
			*Force1 =  (*T.F);					// (u.grad)u-theta			
	}
	
	else if ((globalvar_Pr_switch == "PRSMALL") || (globalvar_Pr_switch == "PRZERO")) {
		
		if (globalvar_RB_Uscaling == "USMALL") 
			*Force1 = (globalvar_Ra)*(*T.F);				// (u.grad)u-theta
		
		else if (globalvar_RB_Uscaling == "ULARGE") 		
			*Force1 = (globalvar_Pr)*(*T.F);				// (u.grad)u-theta
	}
	
	else if (globalvar_Pr_switch == "PRINFTY") {
		*Force1 = (globalvar_Ra)*(*T.F);
	}


	*Force2 = 0.0; 
	*Force3 = 0.0;

	
// For the temperature field

	if (globalvar_Pr_switch == "PRZERO") 
		*T.Force = 0.0;  
	// Do nothing; Nonlinear term (u.grad)T does not exist- see single-time step
	
	else if (globalvar_Pr_switch == "PRLARGE")	
		*T.Force =   globalvar_temperature_grad*(*V1);										
	// F(T) = globalvar_temperature_grad * ux(k)	
	// globalvar_temperature_grad = +1 for RB, and -1 for stratified flows
	
	else if (globalvar_Pr_switch == "PRSMALL") 	
		*T.Force =   complex<DP>(globalvar_temperature_grad/globalvar_Pr, 0)*(*V1); 			
	//F(T) =  globalvar_temperature_grad * ux(k)/Pr
	
	else if (globalvar_Pr_switch == "PRINFTY") {
		*T.Force =   globalvar_temperature_grad*(*V1);
	}


	if (alias_switch == "DEALIAS")		Dealias_force(T);	
	
			
}



//*********************************************************************************************


void IncFluid::Compute_force_RB
(
	IncVF& W, IncSF& T 
)
{
	
	Compute_force_RB(T);

	*W.Force1 = 0.0; 
	*W.Force2 = 0.0; 
	*W.Force3 = 0.0;	

	if (alias_switch == "DEALIAS")		W.Dealias_force();
}


//*********************************************************************************************
//*********************************************************************************************


void IncFluid::Compute_force_RB_rotation(IncSF& T)
{
	// For the velocity field
	
	Compute_force_Coriolis();
	
	// Appropriate multiplication of Coriolis force + buoyancy force
	if (globalvar_Pr_switch == "PRLARGE") {
		
		if (globalvar_RB_Uscaling == "USMALL") 	{
			*Force1 = globalvar_Pr*(*Force1) + (globalvar_Ra*globalvar_Pr)*(*T.F);
			*Force2 = globalvar_Pr*(*Force2);
			*Force3 = globalvar_Pr*(*Force3);
		}
		
		else if (globalvar_RB_Uscaling  == "ULARGE") {
			*Force1 = (sqrt(globalvar_Pr/globalvar_Ra))*(*Force1) + (*T.F);	
			*Force2 = (sqrt(globalvar_Pr/globalvar_Ra))*(*Force2);
			*Force3 = (sqrt(globalvar_Pr/globalvar_Ra))*(*Force3);
		}
	}
	
	else if ((globalvar_Pr_switch == "PRSMALL") || (globalvar_Pr_switch == "PRZERO")) {
		
		if (globalvar_RB_Uscaling == "USMALL") 
			*Force1 = *Force1 + (globalvar_Ra)*(*T.F); // Force2, Force3 unchanged (Coriolos only)		
		
		else if (globalvar_RB_Uscaling == "ULARGE") {
			*Force1 = (sqrt(globalvar_Pr/globalvar_Ra))*(*Force1) + (globalvar_Pr)*(*T.F);	
			*Force2 = (sqrt(globalvar_Pr/globalvar_Ra))*(*Force2);
			*Force3 = (sqrt(globalvar_Pr/globalvar_Ra))*(*Force3);
		}
	}
	
	else if (globalvar_Pr_switch == "PRINFTY") 
		*Force1 = *Force1 + (globalvar_Ra)*(*T.F); // Force2, Force3 unchanged (Coriolos only)		
	
	
	// For the temperature field
	
	if (globalvar_Pr_switch == "PRZERO") 
		*T.Force = 0.0;  
	// Do nothing; Nonlinear term (u.grad)T does not exist- see single-time step
	
	else if (globalvar_Pr_switch == "PRLARGE")	
		*T.Force =   globalvar_temperature_grad*(*V1);										
	// F(T) = globalvar_temperature_grad * ux(k)	
	// globalvar_temperature_grad = +1 for RB, and -1 for stratified flows
	
	else if (globalvar_Pr_switch == "PRSMALL") 	
		*T.Force =   complex<DP>(globalvar_temperature_grad/globalvar_Pr, 0)*(*V1); 			
	//F(T) =  globalvar_temperature_grad * ux(k)/Pr
	
	else if (globalvar_Pr_switch == "PRINFTY") 
	{
		*T.Force =   globalvar_temperature_grad*(*V1);
	}
	
	
	if (alias_switch == "DEALIAS")		Dealias_force(T);	
	
	
}

//*********************************************************************************************
void IncFluid::Compute_force_RB_rotation(IncVF& W, IncSF& T)
{
	Compute_force_RB_rotation(T);
}

//************************ End of compute_force_RB.cc *****************************************

